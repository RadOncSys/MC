#include "mcTransportCylinder.h"
#include "mcGeometry.h"

mcTransportCylinder::mcTransportCylinder()
	:mcTransport()
{
	setGeometry(0, 0);
}

mcTransportCylinder::mcTransportCylinder(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r, double h)
	: mcTransport(orgn, z, x)
{
	setGeometry(r, h);
}

mcTransportCylinder::~mcTransportCylinder(void)
{
}

void mcTransportCylinder::setGeometry(double r, double h)
{
	r_ = r;
	h_ = h;
}

double mcTransportCylinder::getDistanceInside(mcParticle& p) const
{
	return mcGeometry::getDistanceToCylinderInside(p.p, p.u, r_, h_);
}

double mcTransportCylinder::getDistanceOutside(mcParticle& p) const
{
	return mcGeometry::getDistanceToCylinderOutside(p.p, p.u, r_, h_);
}

double mcTransportCylinder::getDNearInside(const geomVector3D& p) const
{
	double cd = r_ - p.lengthXY();
	double dnear = MIN(p.z(), h_ - p.z());
	return MIN(dnear, cd);
}

void mcTransportCylinder::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Radius:\t" << getRadius() << endl;
	os << "Height:\t" << getHeight() << endl;
}

void mcTransportCylinder::dumpVRML(ostream& os) const
{
	dumpVRMLCylinder(os, r_, 0, h_);

	//geomVector3D p = geomVector3D(0, 0, h_*0.5) * mttow_;

	//os << "# Cylinder: " << this->getName() << endl;
	//os << "Transform {" << endl;
	//os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	//os << "  rotation 1 0 0 1.5708" << endl;
	//os << "  children [" << endl;
	//os << "    Shape{" << endl;
	//os << "      appearance Appearance {" << endl;
	//os << "        material Material {" << endl;
	//os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	//os << "          transparency " << transparancy_ << endl;
	//os << "        }" << endl;
	//os << "      }" << endl;
	//os << "      geometry Cylinder { " << endl;
	//os << "                           radius " << r_ << endl;
	//os << "                           height " << h_ << endl;
	//os << "      }" << endl;
	//os << "    }" << endl;
	//os << "  ]" << endl;
	//os << "}" << endl;
}
