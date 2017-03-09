#include "mcSourceSphereC60.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "mcPhysics.h"

mcSourceSphereC60::mcSourceSphereC60(const char* name, int nThreads,
	const geomVector3D& p, const geomVector3D& v, double r)
	:mcSource(name, nThreads)
	, p_(p)
	, v_(v)
	, r_(r)
{
	isGamma_ = true;
}

mcSourceSphereC60::~mcSourceSphereC60(void)
{
}

void mcSourceSphereC60::sample(mcParticle& p, mcThread* thread)
{
	mcRng& rng = thread->rng();

	p.t = MCP_PHOTON;
	p.q = 0;
	p.ke = rng.rnd() < 0.5 ? 1.33 : 1.17;

	// По радиусу распределение пропорционально кубическому корню.
	const static double y = 1.0 / 3.0;
	double r = r_ * pow(rng.rnd(), y);

	mcPhysics::GoInRandomDirection(rng.rnd(), rng.rnd(), p.p);
	p.p *= r;
	p.plast = p.p;

	mcPhysics::GoInRandomDirection(rng.rnd(), rng.rnd(), p.u);

	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke;
}

void mcSourceSphereC60::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	// Источник задается в абсолютных координатах.
	// Поэтому ему не нужна матрица преобразования координат.

	os << "# Source: " << name_ << endl;
	os << "Transform {" << endl;
	os << "  translation " << p_.x() << ' ' << p_.y() << ' ' << p_.z() << endl;
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
