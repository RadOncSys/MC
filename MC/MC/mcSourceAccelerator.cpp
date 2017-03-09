#include "mcSourceAccelerator.h"
#include "mcDefs.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

mcSourceAccelerator::mcSourceAccelerator(const char* name, int nThreads, mc_particle_t type, double ke, double z, double r, double theta)
	:mcSource(name, nThreads)
	, type_(type)
	, q_(0)
	, ke_(ke)
	, z_(z)
	, r_(r)
	, theta_(theta)
{
	tr_ = tan(PI * theta_ / 180);
	uz_ = cos(tr_);
	sinu_ = sin(tr_);
}

void mcSourceAccelerator::sample(mcParticle& p, mcThread* thread)
{
	mcRng& rng = thread->rng();

	p.t = type_;
	p.q = q_;
	p.ke = ke_;

	if (r_ == 0)
		p.p.set(0, 0, z_);
	else
	{
		double phi = rng.rnd() * TWOPI;
		double r = r_ * sqrt(rng.rnd());
		p.p.set(r * cos(phi), r * sin(phi), z_);
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

void mcSourceAccelerator::dumpVRML(ostream& os) const
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
	os << "      geometry Sphere { radius " << r_ << " }" << endl;
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
