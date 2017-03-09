#include "frect.h"
#include <string>

istream& operator >> (istream& is, geomFRect& r)
{
	string line;
	getline(is, line, '\n'); r.l_ = atof(line.c_str());
	getline(is, line, '\n'); r.r_ = atof(line.c_str());
	getline(is, line, '\n'); r.b_ = atof(line.c_str());
	getline(is, line, '\n'); r.t_ = atof(line.c_str());
	return is;
}

ostream& operator << (ostream& os, const geomFRect& r)
{
	os << r.l_ << endl;
	os << r.r_ << endl;
	os << r.b_ << endl;
	os << r.t_ << endl;
	return os;
}
