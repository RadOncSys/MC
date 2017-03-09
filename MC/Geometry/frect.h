// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <math.h>
#include <iostream>
using namespace std;

class geomFRect
{
public:
	geomFRect() :l_(0), r_(0), b_(0), t_(0) {}
	geomFRect(double x, double y)
		:l_(-0.5*x)
		, r_(0.5*x)
		, b_(-0.5*y)
		, t_(0.5*y)
	{}
	geomFRect(double l, double b, double r, double t) :l_(l), r_(r), b_(b), t_(t) {}

	void set(double l, double b, double r, double t) {
		l_ = l; r_ = r; b_ = b; t_ = t;
	}

	inline const double& getL() const { return l_; }
	inline const double& getR() const { return r_; }
	inline const double& getB() const { return b_; }
	inline const double& getT() const { return t_; }

	double	width() const { return fabs(r_ - l_); }
	double	height() const { return fabs(t_ - b_); }
	double	meadleX() const { return (r_ + l_) / 2; }
	double	meadleY() const { return (t_ + b_) / 2; }
	double	minX() const { return (r_ < l_) ? r_ : l_; }
	double	maxX() const { return (r_ > l_) ? r_ : l_; }
	double	minY() const { return (b_ < t_) ? b_ : t_; }
	double	maxY() const { return (b_ > t_) ? b_ : t_; }

	void scale(double f) { l_ *= f; r_ *= f; b_ *= f; t_ *= f; }

	bool isPointInside(double x, double y) const {
		return  x >= minX() && x <= maxX() && y >= minY() && y <= maxY();
	}

	friend istream& operator >> (istream&, geomFRect&);
	friend ostream& operator << (ostream&, const geomFRect&);

public:
	double l_, r_, b_, t_;
};
