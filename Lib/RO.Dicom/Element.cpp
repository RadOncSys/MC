#include "pch.h"
#include "Element.h"
#include "Tag.h"
#include <algorithm>
#include <istream>
#include <locale>
//#include "../../RO.Util/Text.h"

using namespace RO::Dicom;
using namespace std;

const char* unsupported_msg_ = "Dicom::Element: unsuported data type";
const char* wrong_len_msg_ = "element data length do not correspond to element type";
const char* stream_err_msg_ = "error while reading stream";

// Копия старого кода, чтобы развязаться с зависимостью от старого проекта.
short GetTwoStringsFromLine(const string& line, string& s1, string& s2, const char* mask)
{
	s1.erase();
	s2.erase();
	int index = 0, len = (int)line.size();
	if (len == 0) return 0;
	string mStr = (mask == NULL) ? " \t,;=:)(\\" : mask;

	// skip all delimeters in the beginning of line
	while (index < len){
		if (mStr.find_first_of(line[index]) != string::npos)
			index++;
		else
			break;
	}
	if (index >= len) return 0;
	string newline(line, index, len);

	index = (int)newline.find_first_of(mStr);
	if (index == string::npos)   // no sublines (only one line)
		s1 = line;
	else{
		s1 = string(newline, 0, index);
		// skip all delimeters
		len = (int)newline.size();
		while (index < len){
			if (mStr.find_first_of(newline[index]) != string::npos)
				index++;
			else
				break;
		}
		if (index<len)
			s2 = string(newline, index, len);
	}
	return (short)s1.size();
}

int GetFloatArray(const string& line, vector<double>& a)
{
	static const string mask(" \t,;=:)(\\");
	const char* p = line.c_str();
	const char* pend = p + line.size();
	unsigned i = 0;
	bool break_found = true;
	while (true)
	{
		bool is_delimeter = mask.find_first_of(*p) != string::npos;
		if (!break_found || is_delimeter)
		{
			if (is_delimeter) break_found = true;
			p++;
			if (p >= pend) break;
			else continue;
		}
		if (i<a.size()) a[i] = atof(p);
		else a.push_back(atof(p));
		break_found = false;
		p++; i++;
	}
	return i;
}

//
// Utilities
//
UINT16 RO::Dicom::ReadUINT16(std::istream& is)
{
	UINT16 t;
	is.read(reinterpret_cast<char*>(&t), 2);
	return t;
}

UINT32 RO::Dicom::ReadUINT32(std::istream& is)
{
	UINT32 t;
	is.read(reinterpret_cast<char*>(&t), 4);
	return t;
}

UINT32 RO::Dicom::ReadUINT32LE(std::istream& is)
{
	char t[4];
	is.read(t+2, 2);
	is.read(t, 2);
	return *reinterpret_cast<UINT32*>(t);
}

void RO::Dicom::WriteUINT16(std::ostream& os, UINT16 v)
{
	os.write(reinterpret_cast<const char*>(&v), 2);
}

void RO::Dicom::WriteUINT32(std::ostream& os, UINT32 v)
{
	os.write(reinterpret_cast<const char*>(&v), 4);
}

void RO::Dicom::WriteUINT32LE(std::ostream& os, UINT32 v)
{
	const char* p = reinterpret_cast<const char*>(&v);
	os.write(p + 2, 2);
	os.write(p, 2);
}

std::string remove_chars(const std::string& s, const std::string& chars)
{
	std::string result(s.length(), 0);
	auto it = std::copy_if (s.begin(), s.end(), result.begin(),
		[&chars] (const char& c) { return chars.find(c) == std::string::npos; });
	result.resize(std::distance(result.begin(), it));
	return result;
}

std::string remove_non_chars(const std::string& s, const std::string& chars)
{
	std::string result(s.length(), 0);
	auto it = std::copy_if (s.begin(), s.end(), result.begin(),
		[&chars] (const char& c) { return chars.find(c) != std::string::npos; });
	result.resize(std::distance(result.begin(), it));
	return result;
}

std::string remove_non_alnum_chars(const std::string& s, const std::string& chars)
{
	std::locale loc;
	std::string result(s.length(), 0);
	auto it = std::copy_if (s.begin(), s.end(), result.begin(),
		[&chars, loc] (const char& c) { return std::isalnum(c, loc) || chars.find(c) == std::string::npos; });
	result.resize(std::distance(result.begin(), it));
	return result;
}

//
// Element (base class)
//
std::string Element::GetString() const { throw std::exception(unsupported_msg_); }
void Element::SetString(const std::string& s) { throw std::exception(unsupported_msg_); }

std::vector<double> Element::GetDouble() const { throw std::exception(unsupported_msg_); }
void Element::SetDouble(std::vector<double>) { throw std::exception(unsupported_msg_); }

std::vector<int> Element::GetInt() const { throw std::exception(unsupported_msg_); }
void Element::SetInt(std::vector<int>) { throw std::exception(unsupported_msg_); }

std::shared_ptr<DataSet> Element::GetDataSet() const { throw std::exception(unsupported_msg_); }
void Element::SetDataSet(std::shared_ptr<DataSet>) { throw std::exception(unsupported_msg_); }

void Element::GetDate(UINT16& year, UINT16& month, UINT16& day) const { throw std::exception(unsupported_msg_); }
void Element::GetTime(UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const { throw std::exception(unsupported_msg_); }
void Element::GetDateTime(UINT16& year, UINT16& month, UINT16& day, UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const { throw std::exception(unsupported_msg_);}
void Element::SetDate(UINT16 year, UINT16 month, UINT16 day) { throw std::exception(unsupported_msg_); }
void Element::SetTime(UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick) { throw std::exception(unsupported_msg_); }
void Element::SetDateTime(UINT16 year, UINT16 month, UINT16 day, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick) { throw std::exception(unsupported_msg_); }

void Element::ReadTag(std::istream& is) { throw std::exception(unsupported_msg_); }
void Element::ReadVr(std::istream& is) { throw std::exception(unsupported_msg_); }
long Element::ReadLength(std::istream& is) { throw std::exception(unsupported_msg_); }

void Element::WriteTag(std::ostream& os) { throw std::exception(unsupported_msg_); }
void Element::WriteVr(std::ostream& os) { throw std::exception(unsupported_msg_); }
void Element::WriteLength(std::ostream& os, long lengts) { throw std::exception(unsupported_msg_); }

std::wstring Element::ToWString() const
{
	wchar_t txt[64];
	swprintf_s(txt, 64, L"%X\t%X\t%i\t", TagGroup(tag_), TagElement(tag_), (int)vr_);
	return txt;
}

UINT32 Element::PeekUINT32(std::istream& is)
{
	// For some reasons words in the stream are swapped
	char t[4];
	auto pos = is.tellg();
	is.read(t+2, 2);
	is.read(t, 2);
	is.seekg(pos);
	return *reinterpret_cast<UINT32*>(t);
}

UINT32 Element::ReadDelimeter(std::istream& is)
{
	// For some reasons words in the stream are swapped
	char t[4];
	is.read(t+2, 2);
	is.read(t, 2);
	return *reinterpret_cast<UINT32*>(t);
}

//
// String
//
void ElementString::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	data_.resize(len);
	for(UINT32 i=0; i<len; i++)
		data_[i] = is.get();
}

void ElementString::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)data_.size());
	else
		WriteUINT32(os, (UINT16)data_.size());
	os << data_;
}

std::string ElementString::GetString() const { return data_; }
void ElementString::SetString(const std::string& s) { data_ = s; }

std::vector<int> ElementString::GetInt() const
{
	std::vector<int> values;
	string l = data_, l1, l2;
	for (;;){
		if (!GetTwoStringsFromLine(l, l1, l2, nullptr))	break;
		values.push_back(atoi(l1.c_str()));
		l = l2;
	}
	return values;
}

void ElementString::SetInt(std::vector<int> vi)
{
	if (vi.empty()) return;
	char s[32];
	sprintf_s(s, 32, "%i", vi[0]);
	data_ = s;
	for (unsigned i = 1; i < vi.size(); i++)
	{
		data_ += "\\";
		sprintf_s(s, 32, "%i", vi[i]);
		data_ += s;
	}
}

std::vector<double> ElementString::GetDouble() const
{
	std::vector<double> values;
	GetFloatArray(data_, values);
	return values;
}

void ElementString::SetDouble(std::vector<double> vd)
{
	if (vd.empty()) return;
	char s[32];
	sprintf_s(s, 32, "%f", vd[0]);
	data_ = s;
	for (unsigned i = 1; i < vd.size(); i++)
	{
		data_ += "\\";
		sprintf_s(s, 32, "%f", vd[i]);
		data_ += s;
	}
}

std::wstring ElementString::ToWString() const
{
	std::wstring ws(data_.begin(), data_.end());
	return Element::ToWString() + ws;
}

//
// OtherString
//
ElementOtherString::~ElementOtherString()
{
	if(data_ != nullptr)
		delete [] data_;
}

std::string ElementOtherString::GetString() const { return std::string(data_, len_); }
void ElementOtherString::SetString(const std::string& s) { SetData(&s[0], (int)s.length()); }

std::vector<int> ElementOtherString::GetInt() const
{
	std::vector<int> values;
	string l(data_, len_);
	string l1, l2;
	for (;;) {
		if (!GetTwoStringsFromLine(l, l1, l2, nullptr))	break;
		values.push_back(atoi(l1.c_str()));
		l = l2;
	}
	return values;
}

void ElementOtherString::SetInt(std::vector<int> vi)
{
	if (vi.empty()) return;
	char s[32];
	sprintf_s(s, 32, "%i", vi[0]);
	string data = s;
	for (unsigned i = 1; i < vi.size(); i++)
	{
		data += "\\";
		sprintf_s(s, 32, "%i", vi[i]);
		data += s;
	}
	this->SetString(data);
}

std::vector<double> ElementOtherString::GetDouble() const
{
	std::vector<double> values;
	string s(data_, len_);
	GetFloatArray(s, values);
	return values;
}

void ElementOtherString::SetDouble(std::vector<double> vd)
{
	if (vd.empty()) return;
	char s[32];
	sprintf_s(s, 32, "%f", vd[0]);
	string data = s;
	for (unsigned i = 1; i < vd.size(); i++)
	{
		data += "\\";
		sprintf_s(s, 32, "%f", vd[i]);
		data += s;
	}
	this->SetString(data);
}

void ElementOtherString::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	const char* wrong_data = "unexpected delimeter or size";

	// Here argument len, in fact, represents reserved field, which should be 0.
	// So, at this point stream references 32 bit length field.
	// If this length will be equal to 0xffffffff, we have to check delimeter
	// to identify end of element data.

	UINT32 olen = ts->ExplicitVr && len == 0 ? ReadUINT32(is) : len;

	// Undefined length (0xffffffff) can happen inside element for compressed pixel data
	// In this case we should use delimeters
	if(olen == 0xffffffff)
	{
		UINT32 xx = ReadDelimeter(is);
		if(xx != Item) throw std::exception(wrong_data);

		xx = ReadDelimeter(is);
		if(xx != 0) throw std::exception(wrong_data);

		// HACK !!! Догадка. Возможно это уже данные
		xx = ReadDelimeter(is);
		if(xx != Item) throw std::exception(wrong_data);

		// Begin compressed data.
		// The end, probably, defined by
		// This is possible place to put decompressor.

		// Размер здесь будет расширяться динамически
		len_ = 0;
		int maxlen = 512*512*4;
		if(data_ != nullptr)
			delete [] data_;
		data_ = new char[maxlen];

		// Very not efficient.
		// Need better solution by taking the rest of stream.
		// Pixel element should be treated separately from other elements anyway.
		while(true)
		{
			UINT32 nextValue = Element::PeekUINT32(is);
			if(nextValue == SequenceDelimitationItem)
			{
				ReadDelimeter(is);
				xx = ReadDelimeter(is);
				if(xx != 0) throw std::exception(wrong_data);
				break;
			}
			else
			{
				if (len_ > maxlen - 4)
				{
					maxlen += 512*512*4;
					char *ptmp = new char[maxlen];
					memcpy(ptmp, data_, len_);
					delete [] data_;
					data_ = ptmp;
				}
				is.read(data_ + len_, 4);
				if(is.fail()) throw std::exception(stream_err_msg_);
				len_ += 4;
			}
		}
	}
	else
	{
		if(data_ != nullptr)
			delete [] data_;
		data_ = new char[olen];
		is.read(data_, olen);
		len_ = olen;
	}
}

void ElementOtherString::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, 0);
	WriteUINT32(os, len_);
	os.write(data_, len_);
}

void ElementOtherString::SetData(const char* data, int len)
{
	if (data_ != nullptr)
		delete [] data_;
	data_ = new char[len];
	len_ = len;
	memcpy(data_, data, len_);
}

std::wstring ElementOtherString::ToWString() const
{
	std::string s(data_);
	std::wstring ws(s.begin(), s.end());
	return Element::ToWString() + ws;



	//cout << "Group = " << TagGroup(tag) << "      Element = " << TagElement(tag) << endl;



	//if (vd.empty()) return;
	//char s[32];
	//wsprintf_s(s, 32, "%f", vd[0]);
	//string data = s;
	//for (unsigned i = 1; i < vd.size(); i++)
	//{
	//	data += "\\";
	//	sprintf_s(s, 32, "%f", vd[i]);
	//	data += s;
	//}
	//this->SetString(data);




	//auto frmt = safe_cast<Platform::String^>(parameter);
	//auto boxedDouble = dynamic_cast<Box<double>^>(value);
	//auto doubleValue = (boxedDouble == nullptr ? 0 : boxedDouble->Value);
	//wchar_t txt[128];
	//swprintf_s(txt, 128, frmt->Data(), doubleValue);
	//return ref new Platform::String(txt);



	//return std::wstring();


	//return std::wstring();



	//std::string s = e.what();
	//std::wstring ws(s.begin(), s.end());
	//OutputDebugString((ws + L"\n").c_str());
	//throw ref new Platform::FailureException(ref new String(ws.c_str()));
}

//
// AE
//
void ElementAE::SetString(const std::string& s)
{
	// A string of characters that identifies an Application Entity with leading and trailing
	// spaces (20H) being non-significant. A value consisting solely of spaces shall not be used.
	//
	// Default Character Repertoire excluding character code 5CH (the BACKSLASH “\” in ISO-IR 6),
	// and control characters LF, FF, CR and ESC.
	// 16 bytes maximum
	const std::string mask("\\\r\n\f");
	data_ = remove_chars(s, mask);
	if(data_.length() > 16)
		data_.resize(16);
}

//
// AS
//
void ElementAS::SetString(const std::string& s)
{
	const std::string mask("0123456789");
	const std::string mask2("DWMY");
	std::string number = remove_non_chars(s, mask);
	int n = atoi(number.c_str());
	std::string txt = remove_non_chars(s, mask2);
	if(txt.empty())
		throw std::exception("age must contain measurement units symbol");

	char buf[6];
	sprintf_s(buf, 6, "%3u%c", n, txt[0]);
	data_.assign(buf);
}

//
// AT
//
void ElementAT::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	int vm = len / 4;
	for(int i=0; i<vm; i++)
	{
		//UINT16 group = ReadUINT16(is);
		//UINT16 element = ReadUINT16(is);
		//data_.push_back(static_cast<int>(MakeTag(group, element)));
		data_.push_back(static_cast<int>(ReadUINT32LE(is)));
	}
}

void ElementAT::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)(data_.size()*sizeof(int)));
	else
		WriteUINT32(os, (UINT16)(data_.size()*sizeof(int)));
	for (auto&& v : data_)
	{
		WriteUINT32LE(os, v);
	}
}

std::vector<int> ElementAT::GetInt() const { return data_; }
void ElementAT::SetInt(std::vector<int> val) { data_ = val; }

std::wstring ElementAT::ToWString() const
{
	std::wstring ws = Element::ToWString();
	for each (int v in data_)
	{
		wchar_t txt[32];
		swprintf_s(txt, 32, L"%i ", v);
		ws += txt;
	}
	return ws;
}

//
// CS
//
void ElementCS::SetString(const std::string& s)
{
	const std::string mask(" _");
	data_ = remove_non_alnum_chars(s, mask);
	std::transform(data_.begin(), data_.end(), data_.begin(), ::toupper);
	if(data_.length() > 16)
		data_.resize(16);
}

//
// DA
//
void ElementDA::GetDate(UINT16& year, UINT16& month, UINT16& day) const
{
	year = static_cast<UINT16>(atoi(data_.substr(0, 4).c_str()));
	month = static_cast<UINT16>(atoi(data_.substr(4, 2).c_str()));
	day = static_cast<UINT16>(atoi(data_.substr(6, 2).c_str()));
}

void ElementDA::SetDate(UINT16 year, UINT16 month, UINT16 day)
{
	char buf[10];
	sprintf_s(buf, 10, "%04d%02d%02d", year, month, day);
	data_.assign(buf, 8);
}

//
// TM
//
void ElementTM::GetTime(UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const
{
	hour = static_cast<UINT16>(atoi(data_.substr(0, 2).c_str()));
	min = static_cast<UINT16>(atoi(data_.substr(2, 2).c_str()));
	sec = static_cast<UINT16>(atoi(data_.substr(4, 2).c_str()));
	tick = static_cast<UINT16>(atoi(data_.substr(7, 6).c_str()));
}

void ElementTM::SetTime(UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick)
{
	char buf[18];
	sprintf_s(buf, 18, "%04d%02d%02d.%06d ", hour, min, sec, tick);
	data_.assign(buf, 16);
}

//
// DT
//
void ElementDT::GetDateTime(UINT16& year, UINT16& month, UINT16& day, UINT16& hour, UINT16& min, UINT16& sec, UINT32& tick) const
{
	year = static_cast<UINT16>(atoi(data_.substr(0, 4).c_str()));
	month = static_cast<UINT16>(atoi(data_.substr(4, 2).c_str()));
	day = static_cast<UINT16>(atoi(data_.substr(6, 2).c_str()));
	hour = static_cast<UINT16>(atoi(data_.substr(8, 2).c_str()));
	min = static_cast<UINT16>(atoi(data_.substr(10, 2).c_str()));
	sec = static_cast<UINT16>(atoi(data_.substr(12, 2).c_str()));
	tick = static_cast<UINT16>(atoi(data_.substr(15, 6).c_str()));
}

void ElementDT::SetDateTime(UINT16 year, UINT16 month, UINT16 day, UINT16 hour, UINT16 min, UINT16 sec, UINT32 tick)
{
	char buf[28];
	sprintf_s(buf, 10, "%04d%02d%02d%04d%02d%02d.%06d ", year, month, day, hour, min, sec, tick);
	data_.assign(buf, 22);
}

//
// FL
//
void ElementFL::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	int vm = len / 4;
	for(int i=0; i<vm; i++)
	{
		float f;
		is.read(reinterpret_cast<char*>(&f), 4);
		data_.push_back(f);
	}
}

void ElementFL::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)(data_.size()*sizeof(float)));
	else
		WriteUINT32(os, (UINT16)(data_.size()*sizeof(float)));
	for (auto&& v : data_)
	{
		float f = (float) v;
		os.write(reinterpret_cast<char*>(&f), sizeof(float));
	}
}

std::vector<double> ElementFL::GetDouble() const { return data_; }
void ElementFL::SetDouble(std::vector<double> val) { data_ = val; }

std::wstring RO::Dicom::ElementFL::ToWString() const
{
	std::wstring ws = Element::ToWString();
	for each (double v in data_)
	{
		wchar_t txt[32];
		swprintf_s(txt, 32, L"%f ", v);
		ws += txt;
	}
	return ws;
}

//
// FD
//
void ElementFD::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	int vm = len / 8;
	for(int i=0; i<vm; i++)
	{
		double f;
		is.read(reinterpret_cast<char*>(&f), 8);
		data_.push_back(f);
	}
}

void ElementFD::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)(data_.size()*sizeof(double)));
	else
		WriteUINT32(os, (UINT16)(data_.size()*sizeof(double)));
	for (auto v : data_)
	{
		os.write(reinterpret_cast<char*>(&v), sizeof(double));
	}
}

std::vector<double> ElementFD::GetDouble() const { return data_; }
void ElementFD::SetDouble(std::vector<double> val) { data_ = val; }

std::wstring RO::Dicom::ElementFD::ToWString() const
{
	std::wstring ws = Element::ToWString();
	for each (double v in data_)
	{
		wchar_t txt[32];
		swprintf_s(txt, 32, L"%f ", v);
		ws += txt;
	}
	return ws;
}

//
// SL
//
void ElementSL::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	int vm = len / 4;
	for(int i=0; i< vm; i++)
	{
		int f;
		is.read(reinterpret_cast<char*>(&f), 4);
		data_.push_back(f);
	}
}

void ElementSL::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)(data_.size()*sizeof(int)));
	else
		WriteUINT32(os, (UINT16)(data_.size()*sizeof(int)));
	for (auto v : data_)
	{
		os.write(reinterpret_cast<char*>(&v), sizeof(int));
	}
}

std::vector<int> ElementSL::GetInt() const { return data_; }
void ElementSL::SetInt(std::vector<int> val) { data_ = val; }

std::wstring RO::Dicom::ElementSL::ToWString() const
{
	std::wstring ws = Element::ToWString();
	for each (int v in data_)
	{
		wchar_t txt[32];
		swprintf_s(txt, 32, L"%i ", v);
		ws += txt;
	}
	return ws;
}

//
// SS
//
void ElementSS::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	int vm = len / 2;
	for(int i=0; i< vm; i++)
	{
		INT16 f;
		is.read(reinterpret_cast<char*>(&f), 2);
		data_.push_back(static_cast<int>(f));
	}
}

void ElementSS::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)(data_.size()*sizeof(INT16)));
	else
		WriteUINT32(os, (UINT16)(data_.size()*sizeof(INT16)));
	for (auto v : data_)
	{
		os.write(reinterpret_cast<char*>(&v), sizeof(INT16));
	}
}

std::vector<int> ElementSS::GetInt() const { return data_; }
void ElementSS::SetInt(std::vector<int> val) { data_ = val; }

std::wstring RO::Dicom::ElementSS::ToWString() const
{
	std::wstring ws = Element::ToWString();
	for each (int v in data_)
	{
		wchar_t txt[32];
		swprintf_s(txt, 32, L"%i ", v);
		ws += txt;
	}
	return ws;
}

//
// UI
//
void ElementUI::SetString(const std::string& s)
{
	const std::string mask("0123456789.");

	auto it = find_if(s.begin(), s.end(), [mask](const char& c) -> bool
	{
		return mask.find(c) == std::string::npos;
	});

	if( it != s.end() || s.length() > 64)
		throw std::exception("wrong UID format");

	data_.assign(s.c_str());
}

//
// UL
//
void ElementUL::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	int vm = len / 4;
	for(int i=0; i< vm; i++)
	{
		UINT32 f;
		is.read(reinterpret_cast<char*>(&f), 4);
		data_.push_back(static_cast<int>(f));
	}
}

void ElementUL::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)(data_.size()*sizeof(UINT32)));
	else
		WriteUINT32(os, (UINT16)(data_.size()*sizeof(UINT32)));
	for (auto v : data_)
	{
		os.write(reinterpret_cast<char*>(&v), sizeof(UINT32));
	}
}

std::vector<int> ElementUL::GetInt() const { return data_; }
void ElementUL::SetInt(std::vector<int> val) { data_ = val; }

std::wstring RO::Dicom::ElementUL::ToWString() const
{
	std::wstring ws = Element::ToWString();
	for each (int v in data_)
	{
		wchar_t txt[32];
		swprintf_s(txt, 32, L"%i ", v);
		ws += txt;
	}
	return ws;
}

//
// US
//

//template<typename T> const T& RO::Dicom::ElementUS::GetValue() const
//{
//	if(std::is_convertible<UINT16, T>::value)
//		return = dynamic_cast<T>(value);
//	else throw std::exception(unsupported_msg_);
//}
//
//template<typename T> void RO::Dicom::ElementUS::SetValue(T& val)
//{
//	if(std::is_convertible<T, UINT16>::value)
//		data_ = dynamic_cast<UINT16>(value);
//	else throw std::exception(unsupported_msg_);
//}

void ElementUS::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	int vm = len / 2;
	for(int i=0; i< vm; i++)
	{
		UINT16 f;
		is.read(reinterpret_cast<char*>(&f), 2);
		data_.push_back(static_cast<int>(f));
	}
}

void ElementUS::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
		WriteUINT16(os, (UINT16)(data_.size()*sizeof(INT16)));
	else
		WriteUINT32(os, (UINT16)(data_.size()*sizeof(INT16)));
	for (auto v : data_)
	{
		os.write(reinterpret_cast<char*>(&v), sizeof(INT16));
	}
}

std::vector<int> ElementUS::GetInt() const { return data_; }
void ElementUS::SetInt(std::vector<int> val) { data_ = val; }

std::wstring RO::Dicom::ElementUS::ToWString() const
{
	std::wstring ws = Element::ToWString();
	for each (int v in data_)
	{
		wchar_t txt[32];
		swprintf_s(txt, 32, L"%i ", v);
		ws += txt;
	}
	return ws;
}

//
// SQ
//
void ElementSQ::Read(std::istream& is, UINT32 len, TransferSyntaxType ts)
{
	// total sequence length
	unsigned short sz = ReadUINT32(is);		// UINT32 olen
	ReadUINT32(is);		// UINT32 itemDelimeter
	ReadUINT32(is);		// UINT32 itemlen
	if (sz != 0xFFFFFFFF)
		sz -= 8;

	// Read DataSetElements
	data_ = std::make_shared<DataSet>(ts);
	data_->Read(is, sz);
}

void ElementSQ::Write(std::ostream& os, TransferSyntaxType ts) const
{
	if (ts->ExplicitVr)
	{
		WriteUINT16(os, 0);				// reserved
		WriteUINT32LE(os, 0xFFFFFFFF);	// undefined length (sequence)
		WriteUINT32LE(os, Tag::Item);	// item begin delimeter
		WriteUINT32LE(os, 0xFFFFFFFF);	// undefined length (first item)
	}
	else
	{
		WriteUINT32LE(os, 0xFFFFFFFF);
		WriteUINT32LE(os, 0xFFFEE000);
		WriteUINT32LE(os, 0xFFFFFFFF);
		WriteUINT32LE(os, 0x00080000);
	}

	data_->Write(os, ts);

	WriteUINT32LE(os, Tag::ItemDelimitationItem);
	WriteUINT32LE(os, 0x00000000);
	WriteUINT32LE(os, Tag::ItemDelimitationItem);
	WriteUINT32LE(os, 0x00000000);
}

std::shared_ptr<DataSet> ElementSQ::GetDataSet() const { return data_; }
void ElementSQ::SetDataSet(std::shared_ptr<DataSet> ds) { data_ = ds; }
