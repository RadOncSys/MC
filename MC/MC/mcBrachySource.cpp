#include "mcBrachySource.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "mcPhysics.h"

mcBrachySource::mcBrachySource(const char* name, int nThreads,
	const geomVector3D& p, const geomVector3D& v, double r, double h, mc_isotope_t isotope)
	:mcSource(name, nThreads)
	, p_(p)
	, v_(v)
	, r_(r)
	, h_(h)
	, isotope_(isotope)
{
	isGamma_ = true;
}

mcBrachySource::~mcBrachySource(void)
{
}

void mcBrachySource::sample(mcParticle& p, mcThread* thread)
{
	// —пектр Ir-192 представлен количеством фотонов на один распад.
	// —уммарное количество фотонов получаетс€ ё 1, что соответствует физике, 
	// котора€ говорит, что в результате одного распада может по€вл€тьс€ несколько фотонов.
	// ћы здеь не производим перенормировку, что может быть важно в каких-то приложени€х.
	const static double ir192_energies[] = { 0.062245, 0.2058, 0.29596, 0.30847, 0.31651, 0.46807, 0.48458, 0.58859, 0.60442};
	const static double ir192_weights[] = { 0.0523, 0.03303, 0.2867, 0.3269, 0.8286, 0.4783, 0.03187, 0.04515, 0.08232, 0.05309 };

	mcRng& rng = thread->rng();

	p.t = MCP_PHOTON;
	p.q = 0;
	if (isotope_ == mc_isotope_t::IR192)
	{
		int idx = int(rng.rnd() * 10);
		idx = idx < 0 ? 0 : idx > 9 ? 9 : idx;
		p.ke = ir192_energies[idx];
		p.weight = ir192_weights[idx];
	}
	else if (isotope_ == mc_isotope_t::C60)
	{
		p.ke = rng.rnd() < 0.5 ? 1.33 : 1.17;
	}
	else
		throw exception("mcBrachySource; Isotope type not implemeted");

	double z = h_ * rng.rnd();
	double r = r_ * sqrt(rng.rnd());
	double s, c;
	mcPhysics::GetRandomPhi(rng.rnd(), &c, &s);
	p.p.set(p_.x() + r * c, p_.y() + r * s, p_.z() + z);
	p.plast = p.p;

	mcPhysics::GoInRandomDirection(rng.rnd(), rng.rnd(), p.u);

	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke;
}

void mcBrachySource::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	// »сточник задаетс€ в абсолютных координатах.
	// ѕоэтому ему не нужна матрица преобразовани€ координат.

	os << "# Source: " << name_ << endl;
	os << "Transform {" << endl;
	os << "  translation " << p_.x() << ' ' << p_.y() << ' ' << p_.z() + h_ / 2 << endl;
	os << "  rotation 1 0 0 1.5708" << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Cylinder { " << endl;
	os << "                           radius " << r_ << endl;
	os << "                           height " << h_ << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
