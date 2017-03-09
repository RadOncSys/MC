#include "mcTransportSlab.h"
#include <float.h>

mcTransportSlab::mcTransportSlab()
	:mcTransport()
	, h_(0)
{
}

mcTransportSlab::mcTransportSlab(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double h)
	: mcTransport(orgn, z, x)
	, h_(h)
{
}

mcTransportSlab::~mcTransportSlab(void)
{
}

double mcTransportSlab::getDistanceInside(mcParticle& p) const
{
	double vz = p.u.z();
	double pd = (vz < 0) ? -p.p.z() / vz : (vz > 0) ? (h_ - p.p.z()) / vz : DBL_MAX;
	return NNEG(pd);
}

double mcTransportSlab::getDistanceOutside(mcParticle& p) const
{
	// ”даление от секущих плоскостей
	double z = p.p.z(), vz = p.u.z();
	if (z <= 0 && vz <= 0) return DBL_MAX;
	if (z >= h_ && vz >= 0) return DBL_MAX;

	double pd = (vz > 0) ? -p.p.z() / vz : (h_ - p.p.z()) / vz;
	return NNEG(pd);
}

double mcTransportSlab::getDNearInside(const geomVector3D& p) const
{
	double dnear = MIN(p.z(), h_ - p.z());
	return NNEG(dnear);
}

void mcTransportSlab::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Height:\t" << getHeight() << endl;
}

void mcTransportSlab::dumpVRML(ostream& os) const
{
	double a = 50;  // размера бокса
	geomVector3D p = geomVector3D(0, 0, h_*0.5) * mttow_;

	os << "# Slab: " << this->getName() << endl;
	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Box { size " << a << ' ' << a << ' ' << h_ << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
