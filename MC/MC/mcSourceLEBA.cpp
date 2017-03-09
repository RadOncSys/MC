#include "mcSourceLEBA.h"
#include "mcDefs.h"
#include "mcSamplers.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

mcSourceLEBA::mcSourceLEBA(const char* name, int nThreads,
	spectrum_distr_t sptype, profile_distr_t prf_type,
	double kemean, double kewidth, double z, double rsigma)
	:mcSource(name, nThreads)
	, sptype_(sptype)
	, prf_type_(prf_type)
	, q_(-1)
	, kemean_(kemean)
	, kewidth_(kewidth)
	, z_(z)
	, rsigma_(rsigma)
{
}

void mcSourceLEBA::sample(mcParticle& p, mcThread* thread)
{
	mcRng& rng = thread->rng();

	p.t = mc_particle_t::MCP_NEGATRON;
	p.q = q_;

	// Энергия
	double f = 0;
	if (sptype_ == spectrum_distr_t::SPECTRUM_GAUSS)
	{
		// Исключаем из спектра частицы, значения которых за пределами 2 сигма.
		// Для этого линейно переносим диапазон случайных чисел с отрезка (0,1) на (thr, 1-thr).
		//const double thr = 0.005;
		//const double thscale = (1.0 - 2.0 * thr);
		//f = mcSamplers::SampleGauss(thr + rng.rnd() * thscale);
		f = mcSamplers::SampleGauss(rng.rnd());
	}
	else if (sptype_ == spectrum_distr_t::SPECTRUM_PRISM)
	{
		f = mcSamplers::SamplePrism(rng.rnd());
	}
	else if (sptype_ == spectrum_distr_t::SPECTRUM_TRIANGLE)
	{
		f = mcSamplers::SampleTriangle(rng.rnd());
	}
	p.ke = kemean_ + 0.5 * kewidth_ * f;

	// Радиус
	if (rsigma_ == 0)
		p.p.set(0, 0, z_);
	else
	{
		double r = 0;
		if (prf_type_ == profile_distr_t::PROFILE_PRISM)
			r = mcSamplers::SamplePrism(rng.rnd());
		else if (prf_type_ == profile_distr_t::PROFILE_GAUSS)
			r = 2.0 * mcSamplers::SampleGauss(rng.rnd());
		else if (prf_type_ == profile_distr_t::PROFILE_EXPONENT)
			r = mcSamplers::SampleExponent2D(rng.rnd());

		r = fabs(rsigma_ * r);
		double phi = rng.rnd() * TWOPI;
		p.p.set(r * cos(phi), r * sin(phi), z_);
	}
	p.plast = p.p;

	// Направление по определению вдоль оси
	p.u.set(0, 0, 1);
	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke;
}

void mcSourceLEBA::dumpVRML(ostream& os) const
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
	os << "      geometry Sphere { radius " << rsigma_ * sqrt(2.0) << " }" << endl;
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
	os << "                      bottomRadius  " << 0 << endl;
	os << "                      height " << h << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
