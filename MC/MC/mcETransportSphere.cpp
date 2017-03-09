#include "mcETransportSphere.h"
#include "mcGeometry.h"

mcETransportSphere::mcETransportSphere()
	:mcTransport()
	, r_(0)
{
}

mcETransportSphere::mcETransportSphere(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r)
	: mcTransport(orgn, z, x)
	, r_(r)
{
}

double mcETransportSphere::getDistanceInside(mcParticle& p) const
{
	double dist = mcGeometry::getDistanceToSphereInside(p.p, p.u, r_);

	if (internalTransport_ != nullptr)
	{
		// К сожалению частицу нужно переводить в систему координат объекта
		mcParticle pp(p);
		pp.p = pp.p * mtoe_;
		pp.u = p.u.transformDirection(mtoe_);
		double dist2 = internalTransport_->getDistanceOutside(pp);
		if (dist2 < dist)
		{
			p.exitSurface_ = mcParticle::temb_shit_t::Internal;
			return dist2;
		}
		else
		{
			p.exitSurface_ = mcParticle::temb_shit_t::External;
			return dist;
		}
	}
	else
	{
		p.exitSurface_ = mcParticle::temb_shit_t::External;
		return dist;
	}
}

double mcETransportSphere::getDistanceOutside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;
	return mcGeometry::getDistanceToSphereOutside(p.p, p.u, r_);
}

void mcETransportSphere::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Radius:\t" << getRadius() << endl;
}

void mcETransportSphere::dumpVRML(ostream& os) const
{
	mcTransport::dumpVRML(os);

	geomVector3D p = geomVector3D(0, 0, 0) * mttow_;

	os << "# Source: " << name_ << endl;
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
	os << "      geometry Sphere { radius " << r_ << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
