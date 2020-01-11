#include "mcTransportSphereTrap.h"
#include "mcGeometry.h"
#include "mcThread.h"

mcTransportSphereTrap::mcTransportSphereTrap(
	const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r)
	: mcETransportTrap(orgn, z, x), r_(r)
{
}

mcTransportSphereTrap::~mcTransportSphereTrap(void)
{
}

void mcTransportSphereTrap::beginTransport(mcParticle& p)
{
	mcParticle* particle = p.thread_->NextParticle();
	*particle = p;
	particle->p = p.p * mwtot_;
	particle->u = particle->u.transformDirection(mwtot_);
	double R = particle->p.sqLength();
	
	if (R < r_ * r_ && score_)
	{
		double d = mcGeometry::getDistanceToSphereInside(particle->p, particle->u, r_);
		particle->p += particle->u * d;
		score_->ScoreFluence(*particle);
	}
	particle->thread_->RemoveParticle();
}

void mcTransportSphereTrap::dump(ostream& os) const
{
	__super::dump(os);
	os << "Radius:\t" << r_ << endl;
}

void mcTransportSphereTrap::dumpVRML(ostream& os) const
{
	mcTransport::dumpVRML(os);

	geomVector3D p = geomVector3D(0, 0, 0) * mttow_;

	os << "# Transport: " << name_ << endl;
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
