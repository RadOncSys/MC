#include "mcSourceCylindricalC60.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "mcPhysics.h"

mcSourceCylindricalC60::mcSourceCylindricalC60(const char* name, int nThreads,
	const geomVector3D& p, const geomVector3D& v, double r, double h)
	:mcSource(name, nThreads)
	, p_(p)
	, v_(v)
	, r_(r)
	, h_(h)
{
	isGamma_ = true;
}

mcSourceCylindricalC60::~mcSourceCylindricalC60(void)
{
}

void mcSourceCylindricalC60::sample(mcParticle& p, mcThread* thread)
{
	mcRng& rng = thread->rng();

	p.t = MCP_PHOTON;
	p.q = 0;
	p.ke = rng.rnd() < 0.5 ? 1.33 : 1.17;

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

void mcSourceCylindricalC60::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	// Источник задается в абсолютных координатах.
	// Поэтому ему не нужна матрица преобразования координат.

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
