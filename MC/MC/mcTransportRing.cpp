#include "mcTransportRing.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportRing::mcTransportRing(void)
	:mcTransport()
{
	setGeometry(0, 0, 0);
}

mcTransportRing::mcTransportRing(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r0, double r1, double h)
	: mcTransport(orgn, z, x)
{
	setGeometry(r0, r1, h);
}

mcTransportRing::~mcTransportRing(void)
{
}

void mcTransportRing::setGeometry(double r0, double r1, double h)
{
	r0_ = r0;
	r1_ = r1;
	h_ = h;
}

double mcTransportRing::getDistanceInside(mcParticle& p) const
{
	// Цилиндр
	double cd1 = mcGeometry::getDistanceToInfiniteCylinderInside(p.p, p.u, r1_);
	double cd2 = mcGeometry::getDistanceToInfiniteCylinderOutside(p.p, p.u, r0_);
	// Плоскости
	double vz = p.u.z();
	double pd = (vz < 0) ? -p.p.z() / vz : (vz > 0) ? (h_ - p.p.z()) / vz : DBL_MAX;
	return NNEG(MIN(MIN(cd1, cd2), pd));
}

double mcTransportRing::getDistanceOutside(mcParticle& p) const
{
	// Удаление от секущих плоскостей
	double z = p.p.z(), vz = p.u.z();
	if (z <= 0 && vz <= 0) return DBL_MAX;
	if (z >= h_ && vz >= 0) return DBL_MAX;
	double r = p.p.lengthXY();

	// Внутри внутреннего цилиндра
	if (r < r0_)
	{
		double cd = mcGeometry::getDistanceToInfiniteCylinderInside(p.p, p.u, r0_);
		geomVector3D c = p.p + (p.u * cd);
		if (c.z() >= 0 && c.z() <= h_) return cd;
		else if (vz > 0 && c.z() >= h_) return DBL_MAX;
		else if (vz < 0 && c.z() <= 0) return DBL_MAX;
		else if (vz == 0) return DBL_MAX;
		else {
			double pd = (vz > 0) ? -p.p.z() / vz : (h_ - p.p.z()) / vz;
			c = p.p + (p.u * pd);
			r = c.lengthXY();
			if (r > r0_ && r < r1_) return pd;
			else return DBL_MAX;
		}
	}
	// За пределами внешнего цилиндра
	else if (r > r1_)
	{
		double cd = mcGeometry::getDistanceToInfiniteCylinderOutside(p.p, p.u, r1_);
		geomVector3D c = p.p + (p.u * cd);
		if (cd == DBL_MAX) return DBL_MAX;
		else if (c.z() >= 0 && c.z() <= h_) return cd;
		else if (vz > 0 && c.z() >= h_) return DBL_MAX;
		else if (vz < 0 && c.z() <= 0) return DBL_MAX;
		else if (vz == 0) return DBL_MAX;
		else {
			double pd = (vz > 0) ? -p.p.z() / vz : (h_ - p.p.z()) / vz;
			c = p.p + (p.u * pd);
			r = c.lengthXY();
			if (r > r0_ && r < r1_) return pd;
			else return DBL_MAX;
		}
	}
	// Между цилиндрами
	else
	{
		if (vz > 0) {
			double pd = -p.p.z() / vz;
			geomVector3D c = p.p + (p.u * pd);
			r = c.lengthXY();
			if (r > r0_ && r < r1_) return pd;
			else if (r >= r1_) return DBL_MAX;
			else {
				double cd = mcGeometry::getDistanceToInfiniteCylinderInside(c, p.u, r0_);
				c += p.u * cd;
				if (c.z() < h_) return cd + pd;
				else return DBL_MAX; // проскочили в дырку
			}
		}
		else {
			double pd = (h_ - p.p.z()) / vz;
			geomVector3D c = p.p + (p.u * pd);
			r = c.lengthXY();
			if (r > r0_ && r < r1_) return pd;
			else if (r >= r1_) return DBL_MAX;
			else {
				double cd = mcGeometry::getDistanceToInfiniteCylinderInside(c, p.u, r0_);
				c += p.u * cd;
				if (c.z() > 0) return cd + pd;
				else return DBL_MAX;
			}
		}
	}
}

double mcTransportRing::getDNearInside(const geomVector3D& p) const
{
	double z = p.z();
	double r = p.lengthXY();
	double d1 = h_ - z, d2 = r - r0_, d3 = r1_ - r;
	return NNEG(MIN(MIN(z, d1), MIN(d2, d3)));
}

void mcTransportRing::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Radia:\t" << R0() << "\t" << R1() << endl;
	os << "Height:\t" << getHeight() << endl;
}

void mcTransportRing::dumpVRML(ostream& os) const
{
	os << "# Ring: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	dumpVRMLCylinderRing(os, r0_, r1_, 0, h_);

	os << "  ]" << endl;
	os << "}" << endl;
}
