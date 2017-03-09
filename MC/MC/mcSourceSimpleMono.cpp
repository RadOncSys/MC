#include "mcSourceSimpleMono.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

mcSourceSimpleMono::mcSourceSimpleMono(void)
	:mcSource()
	, type_(MCP_PHOTON)
	, ke_(0)
	, q_(0)
{
}

mcSourceSimpleMono::mcSourceSimpleMono(const char* name, int nThreads,
	mc_particle_t type, double ke, const geomVector3D& p, const geomVector3D& v)
	:mcSource(name, nThreads)
{
	init(type, ke, p, v);
}

mcSourceSimpleMono::~mcSourceSimpleMono(void)
{
}

void mcSourceSimpleMono::init(mc_particle_t type
	, double ke
	, const geomVector3D& p
	, const geomVector3D& v)
{
	type_ = type;
	ke_ = ke;
	p_ = p;
	v_ = v;
	q_ = type_ == MCP_PHOTON ? 0 : type_ == MCP_NEGATRON ? -1 : 1;
}

void mcSourceSimpleMono::sample(mcParticle& p, mcThread* thread)
{
	p.t = type_;
	p.q = q_;
	p.ke = ke_;
	p.p = p_;
	p.plast = p.p;
	p.u = v_;
	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += ke_;
}

void mcSourceSimpleMono::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	double r0 = 0.4, r = 0.2, rr = 0.4, h = 1.0, L = 5.0;

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
	os << "      geometry Sphere { radius " << r0 << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;

	geomVector3D p = p_ + (v_ * (L*0.5));

	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
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
	os << "                           radius " << r << endl;
	os << "                           height " << L << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;

	p = p_ + (v_ * (L + h*0.5));

	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
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
	os << "                      bottomRadius  " << rr << endl;
	os << "                      height " << h << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
