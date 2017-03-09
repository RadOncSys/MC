#include "mcETransportTrap.h"
#include "mcScore.h"
#include "mcThread.h"

mcETransportTrap::mcETransportTrap(void)
	:mcTransport()
{
}

mcETransportTrap::mcETransportTrap(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x)
	: mcTransport(orgn, z, x)
{
}

double mcETransportTrap::getDistanceInside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;
	return 0;
}

double mcETransportTrap::getDistanceOutside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;
	return 0;
}

void mcETransportTrap::beginTransport(mcParticle& p)
{
	mcParticle* particle = p.thread_->NextParticle();
	*particle = p;
	particle->p = p.p * mwtot_;
	particle->u = particle->u.transformDirection(mwtot_);
	if (score_)
	{
		score_->ScoreFluence(*particle);
	}
	particle->thread_->RemoveParticle();
}

void mcETransportTrap::beginTransportInside(mcParticle& p)
{
	beginTransport(p);
}

void mcETransportTrap::dumpVRML(ostream& os) const
{
	double a = 30;  // размера бокса
	geomVector3D p = geomVector3D(0, 0, 0.005) * mttow_;

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
	os << "      geometry Box { size " << a << ' ' << a << ' ' << 0.01 << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
