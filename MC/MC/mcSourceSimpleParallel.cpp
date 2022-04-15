#include "mcSourceSimpleParallel.h"
#include "../geometry/vec3d.h"
#include "mcPhysicsCommon.h"
#include "mcThread.h"
#include "mcDefs.h"

mcSourceSimpleParallel::mcSourceSimpleParallel(void)
	:mcSource()
	, type_(MCP_PHOTON)
	, ke_(0)
	, q_(0)
	, rx_(0)
	, ry_(0)
	, sigmaKE_(0)
{
}

mcSourceSimpleParallel::mcSourceSimpleParallel(const char* name, int nThreads,
	mc_particle_t type, double ke, const geomVector3D& p, const geomVector3D& v
	, double rx, double ry, double sigmaKE//=0
)
	:mcSource(name, nThreads)
{
	init(type, ke, p, v, rx, ry, sigmaKE);
}

mcSourceSimpleParallel::~mcSourceSimpleParallel(void)
{
}

void mcSourceSimpleParallel::init(mc_particle_t type
	, double ke
	, const geomVector3D& p
	, const geomVector3D& v
	, double rx				// x-радиус элипса источника
	, double ry				// y-радиус элипса источника
	, double sigmaKE	//= 0.0			// сигма гауссиана энергий
)  // направление движения (единичный вектор)
{
	type_ = type;
	ke_ = ke;
	p_ = p;
	v_ = v;
	q_ = (type_ == mc_particle_t::MCP_NEGATRON) ? -1 : (type_ == mc_particle_t::MCP_POSITRON || type_ == mc_particle_t::MCP_PROTON) ? 1 : 0;
	rx_ = rx;
	ry_ = ry;
	sigmaKE_ = sigmaKE;
}

void mcSourceSimpleParallel::sample(mcParticle& p, mcThread* thread)
{
	p.t = type_;
	p.q = q_;
	mcRng& rng = thread->rng();

	double x, y;
	do {	// единичный круг (x^2+y^2<=1)
		x = 2.0*((rng.rnd()) - 0.5); // положение по x
		y = 2.0*((rng.rnd()) - 0.5); // положение по y
	} while ((SQUARE(x) + SQUARE(y)) > 1);
	x *= rx_; y *= ry_;	// растягиваем в эллипс( (x/rx)^2+(y/ry)^2<=1)
	p.p = geomVector3D(p_.x() + x, p_.y() + y, p_.z()); // не забываем совместить центры координат
	p.plast = p.p;

	if (sigmaKE_ == 0) { p.ke = ke_; }
	else {
		// используем объявленные переменные x,y в другом смысле - для случ. сэмплирования энергии
		// В этом методе сэмплирования ограничение числа сигм видимо ненужно
		GaussStandardRnd_by_Marsaglia(rng, x, y); // случайно распределённые по Гауссу величины
		p.ke = ke_ - sigmaKE_*x;
		p.ke = (p.ke < 0) ? 0 : (p.ke); // проверяем, чтобы энергия не получилась отрицательной
	}

	p.u = v_;
	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += ke_;
}

void mcSourceSimpleParallel::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	os << "# Source: " << name_ << endl;
	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "      transparency " << transparancy_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	os << "        " << p_.x() << ' ' << p_.y() << ' ' << p_.z() << endl;
	os << "        " << p_.x() + v_.x() * 10 << ' ' << p_.y() + v_.y() * 10 << ' ' << p_.z() + v_.z() * 10 << endl;

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	os << "      0 1 -1" << endl;
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;
}
