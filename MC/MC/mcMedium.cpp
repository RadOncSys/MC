#include "mcMedium.h"
#include "../geometry/text.h"

mcMedium::mcMedium(void)
	:status_(EMPTY)
	, density_(0)
{
}

mcMedium::~mcMedium(void)
{
}

string mcMedium::ParseLine(string& line, const char* param)
{
	string s1, s2;
	GetTwoStringsFromLine(line, s1, s2);
	if (s1 != param) throw std::exception((string(param) + string(" expected")).c_str());
	GetTwoStringsFromLine(s2, s1, line);
	return s1;
}
