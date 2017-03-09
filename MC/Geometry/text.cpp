#include "text.h"

using namespace std;

short GetTwoStringsFromLine(const string& line, string& s1, string& s2, const char* mask)
{
	s1.erase();
	s2.erase();
	int index = 0, len = (int)line.size();
	if (len == 0) return 0;
	string mStr = (mask == nullptr) ? " \t,;=:)(\\" : mask;

	// skip all delimeters in the beginning of line
	while (index < len) {
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
	else {
		s1 = string(newline, 0, index);
		// skip all delimeters
		len = (int)newline.size();
		while (index < len) {
			if (mStr.find_first_of(newline[index]) != string::npos)
				index++;
			else
				break;
		}
		if (index < len)
			s2 = string(newline, index, len);
	}
	return (short)s1.size();
}

int GetFloatArray(const string& line, double* x, int n)
{
	static const string mask(" \t,;=:)(\\");
	const char* p = line.c_str();
	const char* pend = p + line.size();
	unsigned i = 0;
	bool break_found = true;
	while (true) {
		bool is_delimeter = mask.find_first_of(*p) != string::npos;
		if (!break_found || is_delimeter) {
			if (is_delimeter) break_found = true;
			p++;
			if (p >= pend) break;
			else continue;
		}
		x[i] = atof(p);
		break_found = false;
		p++; i++;
	}
	return i;
}

int GetFloatArray(const string& line, vector<double>& a)
{
	static const string mask(" \t,;=:)(\\");
	const char* p = line.c_str();
	const char* pend = p + line.size();
	unsigned i = 0;
	bool break_found = true;
	while (true) {
		bool is_delimeter = mask.find_first_of(*p) != string::npos;
		if (!break_found || is_delimeter) {
			if (is_delimeter) break_found = true;
			p++;
			if (p >= pend) break;
			else continue;
		}
		if (i < a.size()) a[i] = atof(p);
		else a.push_back(atof(p));
		break_found = false;
		p++; i++;
	}
	return i;
}

int GetIntArray(const string& line, int* x, int n)
{
	int i = 0;
	string l = line, l1, l2;
	while (i < n && n != 0) {
		if (!GetTwoStringsFromLine(l, l1, l2))
			break;
		x[i] = atoi(l1.c_str());
		l = l2;
		i++;
	}
	return i;
}

int GetIntArray(const string& line, vector<int>& a)
{
	string l = line, l1, l2;
	for (;;) {
		if (!GetTwoStringsFromLine(l, l1, l2))	break;
		a.push_back(atoi(l1.c_str()));
		l = l2;
	}
	return (int)a.size();
}
