#include "mcTransportRectangleRing.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportRectangleRing::mcTransportRectangleRing(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
	double x1, double x2, double y1, double y2, double d, double h)
	:mcTransport(orgn, vz, vx),
	d_(d), h_(h), x1_(x1), x2_(x2), y1_(y1), y2_(y2)
{
}

mcTransportRectangleRing::~mcTransportRectangleRing(void)
{
}

double mcTransportRectangleRing::getDistanceInside(mcParticle& p) const
{
	// Цилиндр
	double cd1 = mcGeometry::getDistanceToRectanglePipeInside(p.p, p.u, x1_ - d_, x2_ + d_, y1_ - d_, y2_ + d_);
	double cd2 = mcGeometry::getDistanceToRectanglePipeOutside(p.p, p.u, x1_, x2_, y1_, y2_);
	// Плоскости
	double vz = p.u.z();
	double pd = (vz < 0) ? -p.p.z() / vz : (vz > 0) ? (h_ - p.p.z()) / vz : DBL_MAX;
	return NNEG(MIN(MIN(cd1, cd2), pd));
}

double mcTransportRectangleRing::getDistanceOutside(mcParticle& p) const
{
	double z = p.p.z(), vz = p.u.z();
	double cd = 0;
	geomVector3D c(p.p);

	// Частица за объектом
	if (z <= 0)
	{
		if (vz <= 0) return DBL_MAX;
		// Доводим частицу до ближайшей плоскости объекта
		cd = -z / vz;
		c = p.p + (p.u * cd);
		if (isXYInside(c.x(), c.y()))
			return cd;
	}

	// Частица перед объектом
	else if (z >= h_)
	{
		if (vz >= 0) return DBL_MAX;
		cd = (h_ - p.p.z()) / vz;
		c = p.p + (p.u * cd);
		if (isXYInside(c.x(), c.y()))
			return cd;
	}

	// Частица на уровне объекта и она либо в дырке, либо за пределами периметра
	double dd = isXYInsideField(c.x(), c.y()) ? mcGeometry::getDistanceToRectanglePipeInside(p.p, p.u, x1_, x2_, y1_, y2_) :
		isXYOutsideCollimator(c.x(), c.y()) ? mcGeometry::getDistanceToRectanglePipeOutside(p.p, p.u, x1_ - d_, x2_ + d_, y1_ - d_, y2_ + d_) : DBL_MAX;

	if (dd == DBL_MAX) return DBL_MAX;
	cd += dd;
	c = p.p + (p.u * cd);
	z = c.z();
	if (z >= 0 && z <= h_) return cd;
	else return DBL_MAX;
}

double mcTransportRectangleRing::getDNearInside(const geomVector3D& p) const
{
	double x = p.x(), y = p.y(), z = p.z();
	double d1 = h_ - z;
	double d2 = distanceToRectangleInside(x, y, x1_ - d_, x2_ + d_, y1_ - d_, y2_ + d_);
	double d3 = distanceToRectangleOutside(x, y, x1_, x2_, y1_, y2_);
	return ZORP(MIN(MIN(z, d1), MIN(d2, d3)));
}

bool mcTransportRectangleRing::isXYInside(double x, double y) const
{
	double X1 = x1_ - d_, X2 = x2_ + d_;
	double Y1 = y1_ - d_, Y2 = y2_ + d_;
	return x > X1 && x < X2 && y > Y1 && y < Y2 &&
		!(x >= x1_ && x <= x2_ && y >= y1_ && y <= y2_);
}

bool mcTransportRectangleRing::isXYInsideField(double x, double y) const
{
	return x >= x1_ && x <= x2_ && y >= y1_ && y <= y2_;
}

bool mcTransportRectangleRing::isXYOutsideCollimator(double x, double y) const
{
	double X1 = x1_ - d_, X2 = x2_ + d_;
	double Y1 = y1_ - d_, Y2 = y2_ + d_;
	return x < X1 || x > X2 || y < Y1 || y > Y2;
}

double mcTransportRectangleRing::distanceToRectangleInside(double x, double y,
	double X1, double X2, double Y1, double Y2) const
{
	return ZORP(MIN(MIN(x - X1, X2 - x), MIN(y - Y1, Y2 - y)));
}

double mcTransportRectangleRing::distanceToRectangleOutside(double x, double y,
	double X1, double X2, double Y1, double Y2) const
{
	if (x <= X1) {
		if (y < Y1) return sqrt((x - X1)*(x - X1) + (y - Y1)*(y - Y1));
		else if (y > Y2) return sqrt((x - X1)*(x - X1) + (y - Y2)*(y - Y2));
		else return ZORP(X1 - x);
	}
	else if (x >= X2) {
		if (y < Y1) return sqrt((x - X2)*(x - X2) + (y - Y1)*(y - Y1));
		else if (y > Y2) return sqrt((x - X2)*(x - X2) + (y - Y2)*(y - Y2));
		else return ZORP(x - X2);
	}
	else {
		if (y < Y1) return ZORP(Y1 - y);
		else if (y > Y2) return ZORP(y - Y2);
		else return 0;
	}
}

void mcTransportRectangleRing::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Height:\t" << h_ << endl;
	os << "Jaw width:\t" << d_ << endl;
	os << "X1 X2:\t" << x1_ << "\t" << x2_ << endl;
	os << "Y1 Y2:\t" << y1_ << "\t" << y2_ << endl;
}

void mcTransportRectangleRing::dumpVRML(ostream& os) const
{
	os << "# Rectangle ring: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	dumpVRMLRectangleRing(os, x1_, x2_, y1_, y2_, d_, h_);

	os << "  ]" << endl;
	os << "}" << endl;
}
