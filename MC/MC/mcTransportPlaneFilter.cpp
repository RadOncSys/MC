#include "mcTransportPlaneFilter.h"
#include "mcScore.h"
#include "mcScoreTrack.h"
#include "mcThread.h"

mcTransportPlaneFilter::mcTransportPlaneFilter(void)
	:mcTransport()
{
}

mcTransportPlaneFilter::mcTransportPlaneFilter(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x)
	: mcTransport(orgn, z, x)
{
}

mcTransportPlaneFilter::~mcTransportPlaneFilter(void)
{
}

void mcTransportPlaneFilter::beginTransport(mcParticle& p)
{
	mcParticle* particle = p.thread_->NextParticle();
	*particle = p;
	particle->p = p.p * mwtot_;
	particle->plast = p.plast * mwtot_;
	particle->u = particle->u.transformDirection(mwtot_);

	if (particle->u.z() > 0 && score_)
	{
		double step = -particle->p.z() / particle->u.z();
		particle->p += particle->u * step;
		score_->ScoreFluence(*particle);
		if (p.trackScore_)
			p.trackScore_->score(p.thread_->id(), p.t, p.p, p.p + (p.u * step), p.ke);
	}

	endTransport(particle);
}

void mcTransportPlaneFilter::dumpVRML(ostream& os) const
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
