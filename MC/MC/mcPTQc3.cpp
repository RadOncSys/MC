#include "mcPTQc3.h"
#include "mcGeometry.h"
#include <float.h>

mcPTQc3::mcPTQc3(void)
	:mcTransport()
{
	setGeometry(0, 0, 0, 0);
}

mcPTQc3::mcPTQc3(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x,
	int nc, double dr, double ds, double h)
	: mcTransport(orgn, z, x)
{
	setGeometry(nc, dr, ds, h);
}

mcPTQc3::~mcPTQc3(void)
{
}

void mcPTQc3::setGeometry(int nc, double dr, double ds, double h)
{
	nc_ = nc;
	dr_ = dr;
	ds_ = ds;
	h_ = h;
}

double mcPTQc3::getDistanceInside(mcParticle& p) const
{
	// Цилиндр
	double r = p.p.lengthXY();
	int ir = int(r / ds_ + MC_EPSILON);
	double r0 = ir * ds_;
	double cd1 = mcGeometry::getDistanceToInfiniteCylinderInside(p.p, p.u, r0 + dr_);
	double cd2 = ir == 0 ? DBL_MAX : mcGeometry::getDistanceToInfiniteCylinderOutside(p.p, p.u, r0);
	// Плоскости
	double vz = p.u.z();
	double pd = (vz < 0) ? -p.p.z() / vz : (vz > 0) ? (h_ - p.p.z()) / vz : DBL_MAX;
	return MIN(MIN(cd1, cd2), fabs(pd));
}

double mcPTQc3::getDistanceOutside(mcParticle& p) const
{
	// Удаление от секущих плоскостей
	double z = p.p.z(), vz = p.u.z();
	if (z <= 0 && vz <= 0) return DBL_MAX;
	if (z >= h_ && vz >= 0) return DBL_MAX;

	// Если частица за пределами торцов, то премещаем ее на ближайщий торец,
	// чтобы определиться с задействованными цилиндрами
	int ir;
	double r, dist = 0;
	geomVector3D pp = p.p;
	if (z <= 0 || z >= h_) {
		dist = (z <= 0) ? fabs(z / vz) : fabs((h_ - z) / vz);
		pp += p.u * dist;
		r = pp.lengthXY();
		ir = int(r / ds_ + MC_EPSILON);
		if (ir < nc_ && r < ir * ds_ + dr_)
			return dist;  // при перенесении на торец воткнулись в полнотелый цылиндр
		z = pp.z();
	}
	else {
		r = pp.lengthXY();
		ir = int(r / ds_ + MC_EPSILON);
	}

	// Определяемся с цилиндрами
	double cd1 = DBL_MAX, cd2 = cd1, cds = cd1;
	if (r > ds_*(nc_ - 1) + dr_)
		cd1 = mcGeometry::getDistanceToInfiniteCylinderOutside(pp, p.u, ds_*(nc_ - 1) + dr_);
	else {
		cd1 = mcGeometry::getDistanceToInfiniteCylinderInside(pp, p.u, (ir + 1) * ds_);
		cd2 = mcGeometry::getDistanceToInfiniteCylinderOutside(pp, p.u, ir * ds_ + dr_);
	}

	// Теперь вопрос, что произойдет быстрее: пересечение с торцом или цилиндром?
	if (vz > 0)
		cds = fabs((h_ - z) / vz);
	else if (vz < 0)
		cds = fabs(z / vz);

	double cd = MIN(cd1, cd2);
	if (cd < cds)
		return dist + cd;
	else
		return DBL_MAX;
}

double mcPTQc3::getDNearInside(const geomVector3D& p) const
{
	double z = p.z();
	double r = p.lengthXY();
	int ir = int(r / ds_ + MC_EPSILON);
	double rr = r - ir * ds_;
	return MIN(MIN(fabs(h_ - z), fabs(z)), MIN(fabs(dr_ - rr), fabs(rr)));
}

void mcPTQc3::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Number of cylinders:\t" << Nc() << endl;
	os << "Radius step:\t" << Ds() << endl;
	os << "Radius thickness:\t" << Dr() << endl;
	os << "Height:\t" << getHeight() << endl;
}

void mcPTQc3::dumpVRML(ostream& os) const
{
	os << "# Cylindrical grid: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	for (int i = 0; i < nc_; i++)
		dumpVRMLCylinderRing(os, i*ds_, i*ds_ + dr_, 0, h_);

	os << "  ]" << endl;
	os << "}" << endl;
}
