#include "mcObj.h"

mcObj::mcObj(void)
	:red_(0.5)
	, green_(0.5)
	, blue_(0.5)
	, transparancy_(0.8)
{
	name_[0] = char(0);
}

void mcObj::setName(const char* s)
{
	strcpy_s(name_, 128, s);
}

void mcObj::setColor(double r, double g, double b, double t)
{
	red_ = r;
	green_ = g;
	blue_ = b;
	transparancy_ = t;
}

ostream& operator << (ostream& os, const mcObj& o)
{
	os << "==========================================================================" << endl;
	os << "Object name:\t" << o.name_ << endl;
	os << "==========================================================================" << endl;
	os << "Color (r,g,b,t):\t" << o.red_ << '\t' << o.green_ << '\t' << o.blue_ << '\t' << o.transparancy_ << endl;
	return os;
}
