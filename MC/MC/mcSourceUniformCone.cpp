#include "mcSourceUniformCone.h"
#include "mcDefs.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

mcSourceUniformCone::mcSourceUniformCone(const char* name, int nThreads, mc_particle_t type, double ke, double z, double r, double sad)
	: mcSource(name, nThreads)
	, type_(type)
	, z_(z)
	, ke_(ke)
	, r_(r)
	, sad_(sad)
{
	transparancy_ = 0.5;
	q_ = (type_ == mc_particle_t::MCP_NEGATRON) ? -1 : (type_ == mc_particle_t::MCP_POSITRON || type_ == mc_particle_t::MCP_PROTON) ? 1 : 0;
}

mcSourceUniformCone::~mcSourceUniformCone(void)
{
}

void mcSourceUniformCone::sample(mcParticle& p, mcThread* thread)
{
	mcRng& rng = thread->rng();
	p.t = type_;
	p.q = q_;
	p.ke = ke_;

	double r = r_ * sqrt(rng.rnd());
	double phi = rng.rnd() * TWOPI;
	double x = r * sin(phi);
	double y = r * cos(phi);

	p.p.set(0, 0, z_);
	p.u.set(x, y, sad_);
	p.u.normalize();

	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke;
}

void mcSourceUniformCone::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	os << "# Uniform cone source: " << this->getName() << endl;
	os << "Transform {" << endl;
	os << "  translation " << "0 " << "0 " << (z_ + sad_ * 0.5) << endl;
	os << "  rotation 1 0 0 " << (-1.5708) << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Cone { " << endl;
	os << "                      bottomRadius " << r_ << endl;
	os << "                      height " << sad_ << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
