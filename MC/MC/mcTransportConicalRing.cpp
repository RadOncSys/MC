#include "mcTransportConicalRing.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportConicalRing::mcTransportConicalRing(void)
	:mcTransport()
{
	setGeometry(0, 0, 0, 0);
}

mcTransportConicalRing::mcTransportConicalRing(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r0, double r1, double h, double f)
	: mcTransport(orgn, z, x)
{
	setGeometry(r0, r1, h, f);
}

mcTransportConicalRing::~mcTransportConicalRing(void)
{
}

void mcTransportConicalRing::setGeometry(double r0, double r1, double h, double f)
{
	r0_ = r0;
	r1_ = r1;
	h_ = h;
	f_ = f;
	if (f > 0)
	{
		cosr0_ = f / sqrt(f*f + r0*r0);
		cosr1_ = f / sqrt(f*f + r1*r1);
	}
	else { cosr0_ = 0; cosr1_ = 0; }
}

double mcTransportConicalRing::getDistanceInside(mcParticle& p) const
{
	// Цилиндр
	double cd1 = mcGeometry::getDistanceToConeInside(p.p, p.u, r1_, f_);
	double cd2 = mcGeometry::getDistanceToConeOutside(p.p, p.u, r0_, f_);
	// Плоскости
	double vz = p.u.z();
	double pd = (vz < 0) ? -p.p.z() / vz : (vz > 0) ? (h_ - p.p.z()) / vz : DBL_MAX;
	return NNEG(MIN(MIN(cd1, cd2), pd));
}

double mcTransportConicalRing::getDistanceOutside(mcParticle& p) const
{
	double z = p.p.z(), vz = p.u.z();
	double cd = 0;
	double rr;	// проекция из фокуса на основание
	geomVector3D c(p.p);

	// Частица за объектом
	if (z <= 0)
	{
		if (vz <= 0) return DBL_MAX;
		// Доводим частицу до ближайшей плоскости объекта
		cd = -p.p.z() / vz;
		c = p.p + (p.u * cd);
		rr = c.lengthXY();
		if (rr < r1_ && rr > r0_) return cd;
	}

	// Частица перед объектом
	else if (z >= h_)
	{
		if (vz >= 0) return DBL_MAX;
		cd = (h_ - p.p.z()) / vz;
		c = p.p + (p.u * cd);
		rr = c.lengthXY() * f_ / (f_ - h_);
		if (rr < r1_ && rr > r0_) return cd;
	}

	else
		rr = c.lengthXY() * f_ / (f_ - z);

	// Частица на уровне объекта и она либо в дырке, либо за пределами кольца
	double dd = rr <= r0_ ? mcGeometry::getDistanceToConeInside(p.p, p.u, r0_, f_) :
		rr >= r1_ ? mcGeometry::getDistanceToConeOutside(p.p, p.u, r1_, f_) : DBL_MAX;

	if (dd == DBL_MAX) return DBL_MAX;
	cd += dd;
	c = p.p + (p.u * cd);
	z = c.z();
	if (z >= 0 && z <= h_) return cd;
	else return DBL_MAX;
}

double mcTransportConicalRing::getDNearInside(const geomVector3D& p) const
{
	double z = p.z();
	double r = p.lengthXY();
	double d1 = h_ - z;
	double d2 = fabs((r - r0_*(f_ - z) / f_) * cosr0_);
	double d3 = fabs((r1_*(f_ - z) / f_ - r) * cosr1_);
	return NNEG(MIN(MIN(z, d1), MIN(d2, d3)));
}

void mcTransportConicalRing::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Radia:\t" << R0() << "\t" << R1() << endl;
	os << "Height:\t" << getHeight() << endl;
	os << "Focus:\t" << F() << endl;
}

void mcTransportConicalRing::dumpVRML(ostream& os) const
{
	os << "# Ring: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	dumpVRMLConicalRing(os, r0_, r1_, 0, h_, f_);

	os << "  ]" << endl;
	os << "}" << endl;
}
