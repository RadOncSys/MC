// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <ole2.h>
#include <string>
#include <vector>

// Аналог узла XML документа.
// В отличие от последнего атрибуты превращаются в узлы, не имеющие детей.
class XPRNode
{
public:
	XPRNode();
	XPRNode(const wchar_t* name, const wchar_t* text);

	XPRNode* AddNode(const wchar_t* name, const wchar_t* text);

	std::wstring Name;
	std::wstring Text;
	std::vector<XPRNode> Nodes;

	// Указатель, необходимый для закрывания узла и возвращения на уровень выше.
	XPRNode* ParentNode;
};

// Класс использует xmllite парсер типа SAX с поддержкой encoding.
// Реализация основана на идее формирования дерева данных по аналогии 
// с DOM, но в сильно упрощенном варианте.
class XmlParseReaderBase
{
public:
	static std::string copyWStringToStlString(const wchar_t* from);
	static inline void ThrowIfFailed(HRESULT hr, const char* s = nullptr);

	static void CreateXPRDocumentFromFile(wchar_t* fname, XPRNode& doc);
	static void CreateXPRDocumentFromStream(IStream * is, XPRNode& doc);

	static std::wstring readElementNodeText(const XPRNode& n, const wchar_t* name);
};
