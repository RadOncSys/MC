// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <iostream>
using namespace std;

class mcObj
{
public:
	mcObj(void);

	void setName(const char*);
	const char* getName() const { return name_; }

	void setColor(double r, double g, double b, double t = 0);
	double r() const { return red_; }
	double g() const { return green_; }
	double b() const { return blue_; }

	friend ostream& operator << (ostream&, const mcObj&);

protected:
	char name_[128];
	double red_;
	double green_;
	double blue_;
	double transparancy_;
};
