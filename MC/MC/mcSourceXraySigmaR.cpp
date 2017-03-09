#include "mcSourceXraySigmaR.h"
#include "mcDefs.h"
#include "mcSamplers.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

mcSourceXraySigmaR::mcSourceXraySigmaR(const char* name, int nThreads, mc_particle_t type, double ke, double z, double sigmar, double theta)
	:mcSource(name, nThreads)
	, type_(type)
	, q_(0)
	, ke_(ke)
	, z_(z)
	, sigma_(sigmar * sqrt(0.5))
	, theta_(theta)
{
	tr_ = tan(PI * theta_ / 180);
	uz_ = cos(tr_);
	sinu_ = sin(tr_);
}

mcSourceXraySigmaR::~mcSourceXraySigmaR(void)
{
}

void mcSourceXraySigmaR::sample(mcParticle& p, mcThread* thread)
{
	// Согласование сигма по радиусу и отдельной координате.
	const double s = sqrt(0.5);
	mcRng& rng = thread->rng();

	p.t = type_;
	p.q = q_;
	p.ke = ke_;

	if (sigma_ == 0)
		p.p.set(0, 0, z_);
	else
	{
		double x = s * sigma_ * mcSamplers::SampleGauss2D(rng.rnd());
		double y = s * sigma_ * mcSamplers::SampleGauss2D(rng.rnd());
		p.p.set(x, y, z_);
	}
	p.plast = p.p;

	if (theta_ == 0)
		p.u.set(0, 0, 1);
	else
	{
		double phi = rng.rnd() * TWOPI;
		p.u.set(sinu_ * cos(phi), sinu_ * sin(phi), uz_);
	}

	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += ke_;
}

void mcSourceXraySigmaR::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	double h = 1.0;

	os << "# Source: " << name_ << endl;

	os << "Transform {" << endl;
	os << "  translation 0 0 " << z_ << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Sphere { radius " << sigma_ * sqrt(2.0) << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;

	os << "Transform {" << endl;
	os << "  translation 0 0 " << z_ - h / 2 << endl;
	os << "  rotation 1 0 0 1.5708" << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Cone { " << endl;
	os << "                      bottomRadius  " << h * sin(tr_) << endl;
	os << "                      height " << h << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
