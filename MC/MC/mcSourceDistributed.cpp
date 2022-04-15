#include "mcSourceDistributed.h"
#include "mcRng.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

extern double const C = 2.99792458e+10;  // cm/s
extern double const EMASS = 510999.06;       // eV
extern double const PMASS = EMASS*1836.1527; // eV
extern double const ECHARGE = 1.6021917e-19;   // Cl
extern double const ALFCURR = 17045.26;        // A
extern double const RECLASS = 2.817939e-13;    // cm

mcSourceDistributed::mcSourceDistributed(void)
	:mcSource(), distr_(), Emitx_(0), rbeamx_(0), rbetax_(0), beamAnglex_(0), Emity_(0), rbeamy_(0), rbetay_(0), beamAngley_(0), beta0_(0)
	, type_(MCP_PHOTON)
	, ke_(0)
	, q_(0)
{
}

mcSourceDistributed::mcSourceDistributed(const char* name, int nThreads,
	mc_particle_t type, double ke, const geomVector3D& p, const geomVector3D& v,
	mc_distr_t distr, double Emitx, double rbeamx, double beamAnglex, double Emity,
	double rbeamy, double beamAngley)
	:mcSource(name, nThreads)
{
	init(type, ke, p, v, distr, Emitx, rbeamx, beamAnglex, Emity, rbeamy, beamAngley);
}

void mcSourceDistributed::init(mc_particle_t type
	, double ke
	, const geomVector3D& p
	, const geomVector3D& v
	, mc_distr_t distr        // тип распределения
	, double Emitx  // х-эмиттанс, в см. "Нормализованный", без пи, т.е. просто rbeam*rbeta, rbeta = vrmax/C, C - скорость света
	, double rbeamx // полуось по х
	, double beamAnglex // угловое расхождение по х (из-за оптики, не связано с эмиттансом)
	, double Emity            //
	, double rbeamy           //  все то же самое для y
	, double beamAngley       //
)
{
	type_ = type;
	ke_ = ke;
	p_ = p;
	v_ = v;
	q_ = (type_ == mc_particle_t::MCP_NEGATRON) ? -1 : (type_ == mc_particle_t::MCP_POSITRON || type_ == mc_particle_t::MCP_PROTON) ? 1 : 0;

	distr_ = distr;  // тип распределения
	Emitx_ = Emitx;  // х-эмиттанс, в см. "Нормализованный", без пи, т.е. просто rbeam*rbeta, rbeta = vrmax/C, C - скорость света
	rbeamx_ = rbeamx; // полуось по х
	beamAnglex_ = beamAnglex; // угловое расхождение по х (из-за оптики, не связано с эмиттансом)
	Emity_ = Emity;                   //
	rbeamy_ = rbeamy;                //  все то же самое для y
	beamAngley_ = beamAngley;       //

	beta0_ = sqrt(ke_*1e6 / PMASS / (ke_*1e6 / PMASS + 1)); // если источник не протонный - исправить массу!!!!! Предполагается, что ke_ задано в МэВ-ах
	rbetax_ = Emitx_ / rbeamx_;
	rbetay_ = Emity_ / rbeamy_;

	//  fdistr_ = NULL;  // для произвольного распределения дать указатель на внешнюю функцию распределения.
}

void mcSourceDistributed::sample(mcParticle& p, mcThread* thread)
{
	double x = 0, y = 0, px = 0, py = 0;  //координаты и импульсы
	mcRng& rng = thread->rng();

	// если однородный по плотности или гауссов - готовим однородный круг в x,y
	if (distr_ == MCP_CONST_DENS || distr_ == MCP_GAUSSIAN) {
		double r2, rbeta2, spread;
		do {
			x = (rng.rnd() - 0.5)*rbeamx_ * 2;
			y = (rng.rnd() - 0.5)*rbeamy_ * 2;
		} while (x*x / rbeamx_ / rbeamx_ + y*y / rbeamy_ / rbeamy_ > 1);

		// если гауссов  - растягиваем созданный круг

		static double width = 3.;   // на скольких сигмах обрезка 
		static double mult = sqrt(1 - exp(-width*width)); // сначала немного однородно сжимаем
		if (distr_ == MCP_GAUSSIAN) {
			x *= mult;
			y *= mult;
			r2 = x*x / rbeamx_ / rbeamx_ + y*y / rbeamy_ / rbeamy_;
			spread = sqrt(-log(1 - r2) / r2); // и потом неоднородно растягиваем
			x *= spread;
			y *= spread;
		};

		//теперь все то же самое с импульсами, но с учетом квадратичного уменьшения температуры от центра к границе

		//для однородной плотности
		do {
			px = (rng.rnd() - 0.5) * 2 * rbetax_;
			py = (rng.rnd() - 0.5) * 2 * rbetay_;
		} while (px*px / rbetax_ / rbetax_ + py*py / rbetay_ / rbetay_
				> (1 - x*x / rbeamx_ / rbeamx_ - y*y / rbeamy_ / rbeamy_));

		//для гауссова пучка
		if (distr_ == MCP_GAUSSIAN) {
			do {
				px = (rng.rnd() - 0.5) * 2 * rbetax_;
				py = (rng.rnd() - 0.5) * 2 * rbetay_;
			} while (px*px / rbetax_ / rbetax_ + py*py / rbetay_ / rbetay_
				  > (1 - x*x / rbeamx_ / rbeamx_ / width / width - y*y / rbeamy_ / rbeamy_ / width / width));
			px *= mult;
			py *= mult;
			rbeta2 = px*px / rbetax_ / rbetax_ + py*py / rbetay_ / rbetay_;
			spread = sqrt(-log(1 - rbeta2) / rbeta2);
			px *= spread;
			py *= spread;
		};
	};

	// если однородный 4-мерный эллипсоид (он же waterbag, "водяной мешок")
	if (distr_ == MCP_WATERBAG || distr_ == MCP_K_V) {
		double Rx2, Ry2;
		do {
			do {
				x = (rng.rnd() - 0.5)*rbeamx_ * 2;
				px = (rng.rnd() - 0.5) * 2 * rbetax_;
			} while ((Rx2 = x*x / rbetax_ / rbetax_ + px*px / rbetax_ / rbetax_) > 1);
			do {
				y = (rng.rnd() - 0.5)*rbeamy_ * 2;
				py = (rng.rnd() - 0.5) * 2 * rbetay_;
			} while ((Ry2 = y*y / rbetay_ / rbetay_ + py*py / rbetay_ / rbetay_) > 1);
		} while ((Rx2 + Ry2) > 1);

		// waterbag легко преобразуется к распределению Капчинского-Владимирского
		// т.е. к однородной тонкой оболочке 4-мерного эллипсоида.
		// В проекциях на фазовые плоскости (x,px) и (y,py) получаются эллипсы с однородной плотностью,
		// довольно часто применяется при моделировании
		if (distr_ == MCP_K_V) {
			double rr = sqrt(Rx2 + Ry2);
			if (rr) {
				x /= rr;
				px /= rr;
				y /= rr;
				py /= rr;
			};
		};
	};

	if (distr_ == MCP_ARBITRARY) { //  произвольное распределение, пока только демонстрация алшоритма
		double Rx2, Ry2;
		do {
			do {
				x = (rng.rnd() - 0.5)*rbeamx_ * 2;
				px = (rng.rnd() - 0.5) * 2 * rbetax_;
			} while ((Rx2 = x*x / rbetax_ / rbetax_ + px*px / rbetax_ / rbetax_) > 1);
			do {
				y = (rng.rnd() - 0.5)*rbeamy_ * 2;
				py = (rng.rnd() - 0.5) * 2 * rbetay_;
			} while ((Ry2 = y*y / rbetay_ / rbetay_ + py*py / rbetay_ / rbetay_) > 1);
			//для работы - убрать комментарий и в init определить ссылку на внешнюю функцию *fdistr_
		} while ((Rx2 + Ry2) > 1 /* && rng.rnd() > fdistr_(x,y,px,py) */);
	};

	// напоследок - учет углового расхождения пучка
	px += beta0_*beamAnglex_ / rbeamx_*x;
	py += beta0_*beamAngley_ / rbeamy_*y;

	// готовим частицу
	p.t = type_;
	p.q = q_;
	p.ke = ke_;
	p.p = p_ + geomVector3D(x, y, 0);
	// углы считаем малыми
	p.u = geomVector3D(px / beta0_, py / beta0_, sqrt(1. - (px / beta0_*px / beta0_) - (py / beta0_*py / beta0_)));
	p.weight = 1;     // для ионов - изменить
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += ke_;
}

void mcSourceDistributed::dumpVRML(ostream& os) const
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
