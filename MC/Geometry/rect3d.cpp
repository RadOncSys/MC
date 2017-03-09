#include "rect3d.h"

geomRect3D::geomRect3D()
{
}

geomRect3D::geomRect3D(const geomRect3D& r)
{
	set(r.getPoint()
		, r.getNormal()
		, r.getXAxis()
		, (const geomFRect&)r);
}

void geomRect3D::set(const geomVector3D& p
	, const geomVector3D& n
	, const geomVector3D& xv
	, const geomFRect& rect)
{
	geomPlane3D::set(p, n, xv);
	*(geomFRect*)this = rect;
}

geomVector3D geomRect3D::getCenter()const
{
	return p_ + geomVector3D(meadleX(), meadleY(), 0);
}

void geomRect3D::operator=(const geomRect3D& r)
{
	set(r.getPoint()
		, r.getNormal()
		, r.getXAxis()
		, (const geomFRect&)r);
}

istream& operator >> (istream& is, geomRect3D& r)
{
	is >> (geomPlane3D&)r;
	is >> (geomFRect&)r;
	return is;
}

ostream& operator << (ostream& os, const geomRect3D& r)
{
	os << (geomPlane3D&)r;
	os << (geomFRect&)r;
	return os;
}
