#include "plane3d.h"

geomPlane3D::
geomPlane3D()
{
}

geomPlane3D::
geomPlane3D(const geomPlane3D& p)
	:p_(p.getPoint())
	, n_(p.getNormal())
	, xv_(p.getXAxis())
	, mRtoP_(p.getRtoPMatrix())
	, mPtoR_(p.getPtoRMatrix())
{}

void
geomPlane3D::
set(const geomVector3D& p, const geomVector3D& n, const geomVector3D& xv)
{
	/*
	// HACK !!!
	// Из-за ошибок округления возникают проблемы сопоставления дозовых матриц со срезами.
	// Здесь точка плоскости окргляется до 3-го знака (эквивалентно 0.01 мм).
	// В других приложениях тип микромира :) или очень больших чисел (> INT_MAX/1000) могут быть проблемы.
	p_.set(0.001*ROUND(p.x()*1000.)
	,0.001*ROUND(p.y()*1000.)
	,0.001*ROUND(p.z()*1000.));
	*/
	p_ = p;
	n_ = n;
	xv_ = xv;
	n_.normalize();
	xv_.normalize();
	initMRtoP();
}

void
geomPlane3D::
initMRtoP()
{
	geomVector3D yv = n_ ^ xv_;		// y - axis of the plane
	static const geomVector3D i(1, 0, 0);
	static const geomVector3D j(0, 1, 0);
	static const geomVector3D k(0, 0, 1);

	// Beam to Patient Translation
	mRtoP_(0, 0) = i * xv_; mRtoP_(0, 1) = j * xv_;	mRtoP_(0, 2) = k * xv_; mRtoP_(0, 3) = 0;
	mRtoP_(1, 0) = i * yv;	 mRtoP_(1, 1) = j * yv;	mRtoP_(1, 2) = k * yv;	 mRtoP_(1, 3) = 0;
	mRtoP_(2, 0) = i * n_;	 mRtoP_(2, 1) = j * n_;	mRtoP_(2, 2) = k * n_;	 mRtoP_(2, 3) = 0;
	mRtoP_(3, 0) = p_(0);	 mRtoP_(3, 1) = p_(1);	  mRtoP_(3, 2) = p_(2);	 mRtoP_(3, 3) = 1;

	mPtoR_ = mRtoP_;
	mPtoR_.makeInverse();
}

geomVector3D
geomPlane3D::
crossByLine(const geomVector3D& p0, const geomVector3D& p1)const
{
	double h0 = (p0 - p_)*n_, h1 = (p1 - p_)*n_;
	if (h0 == h1)
		throw std::exception("geomPlane3D:: линия параллельна плоскости");
	return p0 + (p1 - p0)*(h0 / (h0 - h1));
}

double
geomPlane3D::
getPlanePosition() const
{
	return p_ * n_;
}

bool
geomPlane3D::
crossByEdge(const geomVector3D& p0, const geomVector3D& p1, double& x, double& y)const
{
	double h0 = (p0 - p_)*n_, h1 = (p1 - p_)*n_;
	if (h0*h1 > 0) return false;
	double dh = h0 - h1;
	geomVector3D pc = dh == 0 ? p0 : p0 + (p1 - p0)*(h0 / dh);
	pc = pc * mPtoR_;
	x = pc.x();
	y = pc.y();
	return true;
}

double
geomPlane3D::
nearestDistance(const geomVector3D& p)const
{
	return fabs((p - p_)*n_);
}

istream& operator >> (istream& is, geomPlane3D& plane)
{
	is >> plane.p_;
	is >> plane.n_;
	is >> plane.xv_;
	plane.initMRtoP();
	return is;
}

ostream& operator << (ostream& os, const geomPlane3D& plane)
{
	os << plane.p_;
	os << plane.n_;
	os << plane.xv_;
	return os;
}
