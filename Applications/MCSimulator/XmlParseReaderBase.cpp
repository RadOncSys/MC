#include "XmlParseReaderBase.h"
#include <xmllite.h>
#include <comdef.h>
#include <shlwapi.h>

using namespace std;

#define SAFE_RELEASE(I)         do { if (I){ I->Release(); } I = NULL; } while(0)

inline void XmlParseReaderBase::ThrowIfFailed(HRESULT hr, const char* s)
{
	if (FAILED(hr))
	{
		if (s == nullptr)
		{
			auto msg = _com_error(hr).ErrorMessage();
			auto len = lstrlenW(msg);
			std::string str((size_t)len + 1, 0);
			WideCharToMultiByte(CP_ACP, 0, msg, -1, &str[0], len, nullptr, nullptr);
			throw std::exception(str.c_str());
		}
		else
			throw std::exception(s);
	}
}

XPRNode::XPRNode() : ParentNode(nullptr) {}

XPRNode::XPRNode(const wchar_t* name, const wchar_t* text) : ParentNode(nullptr)
{
	if (name != nullptr) Name = name;
	if (text != nullptr) Text = text;
}

XPRNode* XPRNode::AddNode(const wchar_t* name, const wchar_t* text)
{
	XPRNode node(name, text);
	node.ParentNode = this;
	Nodes.push_back(node);
	return &Nodes.back();
}

string XmlParseReaderBase::copyWStringToStlString(const wchar_t* from)
{
	if (from == nullptr) return string();
	auto len = lstrlenW(from);
	std::string str((size_t)len + 1, 0);
	WideCharToMultiByte(CP_ACP, 0, from, -1, &str[0], len, nullptr, nullptr);
	return str;
}

void XmlParseReaderBase::CreateXPRDocumentFromFile(wchar_t* fname, XPRNode& doc)
{
	string errmsg = "Error creating XML reader for file: " + copyWStringToStlString(fname);
	IStream *pFileStream = nullptr;
	ThrowIfFailed(SHCreateStreamOnFile(fname, STGM_READ, &pFileStream), errmsg.c_str());
	CreateXPRDocumentFromStream(pFileStream, doc);
	SAFE_RELEASE(pFileStream);
}

void XmlParseReaderBase::CreateXPRDocumentFromStream(IStream* is, XPRNode& doc)
{
	// Создание COM обекта - ридер XML
	IXmlReader *reader = nullptr;
	XmlParseReaderBase::ThrowIfFailed(CreateXmlReader(
		__uuidof(IXmlReader), (void**)&reader, nullptr), "Error creating xml reader");

	// Режима работы ридера
	XmlParseReaderBase::ThrowIfFailed(
		reader->SetProperty(XmlReaderProperty_DtdProcessing, DtdProcessing_Prohibit), "Error setting XmlReaderProperty_DtdProcessing");

	// Подключение к потоку данных
	XmlParseReaderBase::ThrowIfFailed(reader->SetInput(is), "Error setting input for reader");

	XmlNodeType nodeType;
	LPCWSTR pwszPrefix = nullptr;
	LPCWSTR pwszLocalName = nullptr;
	LPCWSTR pwszValue = nullptr;
	UINT cwchPrefix = 0, cwchLocalName = 0, cwchValue = 0;
	XPRNode* currentNode = nullptr;
	UINT currentDepth = 0, nodeDepth = 0, attrCount = 0;

	// Reads until there are no more nodes.
	while (S_OK == reader->Read(&nodeType))
	{
		switch (nodeType)
		{
		case XmlNodeType_XmlDeclaration:
			//text += L"XmlDeclaration\n";
			//ThrowIfFailed(ReadAttributes(reader.Get(), text));
			break;

		case XmlNodeType_Element:
			ThrowIfFailed(reader->GetPrefix(&pwszPrefix, &cwchPrefix));
			ThrowIfFailed(reader->GetLocalName(&pwszLocalName, &cwchLocalName));
			ThrowIfFailed(reader->GetDepth(&nodeDepth));
			ThrowIfFailed(reader->GetAttributeCount(&attrCount));

			// Если это корневой узел, то просто читаем в него.
			// в противном случае создаем новый вложенный.
			if (nodeDepth == 0)
			{
				currentNode = &doc;
				currentDepth = 0;
			}
			else if (nodeDepth > currentDepth)
			{
				currentNode = currentNode->AddNode(nullptr, nullptr);
				currentDepth = nodeDepth;
			}
			else
			{
				currentNode = currentNode->ParentNode->AddNode(nullptr, nullptr);
				currentDepth = nodeDepth;
			}
			currentNode->Name = wstring(pwszPrefix, 0, cwchPrefix) + wstring(pwszLocalName, 0, cwchLocalName);

			// Read attributes
			if (attrCount > 0)
			{
				ThrowIfFailed(reader->MoveToFirstAttribute());
				while (true)
				{
					ThrowIfFailed(reader->GetPrefix(&pwszPrefix, &cwchPrefix));
					ThrowIfFailed(reader->GetLocalName(&pwszLocalName, &cwchLocalName));
					ThrowIfFailed(reader->GetValue(&pwszValue, &cwchValue));

					if (cwchValue > 0)
					{
						currentNode->AddNode(
							(wstring(pwszPrefix, 0, cwchPrefix) + wstring(pwszLocalName, 0, cwchLocalName)).c_str(),
							wstring(pwszValue, 0, cwchValue).c_str());
					}
					if (S_OK != reader->MoveToNextAttribute())
						break;
				}
			}
			break;

		case XmlNodeType_EndElement:
			ThrowIfFailed(reader->GetPrefix(&pwszPrefix, &cwchPrefix));
			ThrowIfFailed(reader->GetLocalName(&pwszLocalName, &cwchLocalName));
			currentNode = currentNode->ParentNode;
			currentDepth--;
			break;

		case XmlNodeType_Whitespace:
			//ThrowIfFailed(reader->GetValue(&pwszValue, &cwchValue));
			//text += L"Whitespace: >";
			//text += ref new String(pwszValue, cwchValue);
			//text += L"<\n";
			break;

		case XmlNodeType_Text:
			ThrowIfFailed(reader->GetValue(&pwszValue, &cwchValue));
			currentNode->Text = wstring(pwszValue, 0, cwchValue);
			break;

		case XmlNodeType_CDATA:
			//ThrowIfFailed(reader->GetValue(&pwszValue, &cwchValue));
			//text += L"CDATA: ";
			//text += ref new String(pwszValue, cwchValue);
			//text += L"\n";
			break;
		case XmlNodeType_ProcessingInstruction:
			//ThrowIfFailed(reader->GetLocalName(&pwszLocalName, &cwchLocalName));
			//ThrowIfFailed(reader->GetValue(&pwszValue, &cwchValue));
			//text += L"Processing Instruction name:";
			//text += ref new String(pwszLocalName, cwchLocalName);
			//text += L"value:";
			//text += ref new String(pwszValue, cwchValue);
			//text += L"\n";
			break;
		case XmlNodeType_Comment:
			//ThrowIfFailed(reader->GetValue(&pwszValue, &cwchValue));
			//text += L"Comment: ";
			//text += ref new String(pwszValue, cwchValue);
			//text += L"\n";
			break;
		case XmlNodeType_DocumentType:
			//text += L"DOCTYPE is not printed\n";
			break;
		case XmlNodeType_None: break;
		case XmlNodeType_Attribute: break;
		default: break;
		}
	}
	SAFE_RELEASE(reader);
}

wstring XmlParseReaderBase::readElementNodeText(const XPRNode& n, const wchar_t* name)
{
	wstring text;
	for (auto node : n.Nodes)
	{
		if (node.Name == name)
		{
			text = node.Text;
			break;
		}
	}
	return text;
}
