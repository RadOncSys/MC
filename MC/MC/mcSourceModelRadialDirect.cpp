#include "mcSourceModelRadialDirect.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "mcPhysics.h"
#include "mcDefs.h"
#include <math.h>

mcSourceModelRadialDirect::mcSourceModelRadialDirect(const char* name, int nThreads, double z0)
	:mcSource(name, nThreads), z0_(z0), particles_(nullptr)
{
	currentParticleIdx_.resize(nThreads_);
	SetSplitting(1);
}

mcSourceModelRadialDirect::~mcSourceModelRadialDirect()
{
	if (particles_ != nullptr)
		delete [] particles_;
}

void mcSourceModelRadialDirect::sample(mcParticle& p, mcThread* thread)
{
	auto rng = thread->rng();
	auto tidx = thread->id();

	p.t = MCP_PHOTON;
	p.q = 0;

	auto pidx = currentParticleIdx_[tidx];
	auto sidx = currentSplitIdx_[tidx];
	currentSplitIdx_[tidx]++;

	// Переходим к новой частице по окнчании пула расщепленных текущей.
	if (noffsplits_ == 0 || currentSplitIdx_[tidx] % noffsplits_ == 0)
	{
		currentSplitIdx_[tidx] = 0;
		currentParticleIdx_[tidx]++;
		if (currentParticleIdx_[tidx] % noofParticles_ == 0)
			currentParticleIdx_[tidx] = 0;
		splitStartAngle_[tidx] = TWOPI * rng.rnd();
	}

	const unsigned short* sp = particles_ + pidx * 4;

	p.plast = p.p;
	p.weight = dw_;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;

	p.ke = sp[0] * energyScale_;

	double r = sp[1] * rScale_;
	double ux = sp[2] * radialAngleScale_ - 1.0, uy = sp[3] * azimutAngleScale_ - 1.0, uz = sqrt(1.0 - ux * ux - uy * uy);

	// Поворачиваем вокруг оси системы.
	double phi = splitStartAngle_[tidx] + sidx * da_;
	double cosPhi = cos(phi);
	double sinPhi = sin(phi);

	p.p.set(r, 0, 0);
	p.p.p_[0] = r * cosPhi;
	p.p.p_[1] = -r * sinPhi;
	p.p.p_[2] += z0_;

	p.u.p_[0] = ux * cosPhi + uy * sinPhi;
	p.u.p_[1] = -ux * sinPhi + uy * cosPhi;
	p.u.p_[2] = uz;

	etotal_[tidx] += p.ke * p.weight;
}

void mcSourceModelRadialDirect::Init(double maxEnergy, double maxR)
{
	energyScale_ = maxEnergy / 0xffff;
	rScale_ = maxR / 0xffff;
	radialAngleScale_ = 2.0 / 0xffff;
	azimutAngleScale_ = 2.0 / 0xffff;
}

void mcSourceModelRadialDirect::SetSplitting(unsigned nSplits)
{
	noffsplits_ = nSplits;
	da_ = nSplits == 0 ? 0 : TWOPI / nSplits;
	dw_ =  nSplits == 0 ? 1.0 : 1.0 / nSplits;
	splitStartAngle_.resize(nThreads_, 0);
	currentSplitIdx_.resize(nThreads_, 0);
}

void mcSourceModelRadialDirect::SetParticles(unsigned short* data, unsigned nThreadSize, const std::vector<unsigned>& indexes)
{
	noofParticles_ = 0;
	for (unsigned i = 0; i < indexes.size(); i++)
		noofParticles_ += indexes[i];

	particles_ = new unsigned short[noofParticles_ * 4];
	unsigned short* pdest = particles_;

	unsigned idx = 0;
	for (unsigned i = 0; i < indexes.size(); i++)
	{
		const unsigned short* psrc = data + nThreadSize * i * 4;
		for (unsigned j = 0; j < indexes[i] * 4; j++, psrc++, pdest++)
			*pdest = *psrc;
	}
}

void* mcSourceModelRadialDirect::saveToMemory(int& size)
{
	int nparsBytes = 5 * sizeof(double) + 1 * sizeof(unsigned);
	int nparticlesBytes = noofParticles_ * 4 * sizeof(unsigned short);
	size = nparsBytes + nparticlesBytes;

	void* buffer = malloc(size);
	if (buffer == nullptr)
		throw std::exception("mcSourceModelRadialPhsp::saveToMemory can not allocate memory");

	char* pbuffer = (char*)buffer;
	*(double*)pbuffer = z0_; pbuffer += sizeof(double);
	*(double*)pbuffer = energyScale_; pbuffer += sizeof(double);
	*(double*)pbuffer = rScale_; pbuffer += sizeof(double);
	*(double*)pbuffer = radialAngleScale_; pbuffer += sizeof(double);
	*(double*)pbuffer = azimutAngleScale_; pbuffer += sizeof(double);
	*(unsigned*)pbuffer = noofParticles_; pbuffer += sizeof(unsigned);

	memcpy(pbuffer, particles_, nparticlesBytes);

	return buffer;
}

void mcSourceModelRadialDirect::readFromMemory(void* buffer)
{
	char* pbuffer = (char*)buffer;

	// HACK!! Несколько грязно.
	// Положение источника есть в модели и в конструкторю.
	// Пока значение в файле модели игнорируем.
	// Используем значение в конструкторе.
	// Это опасный момент и надо с ним разобраться позднее.
	double z0;

	// z0_ = *(double*)(pbuffer); pbuffer += sizeof(double);
	z0 = *(double*)(pbuffer); pbuffer += sizeof(double);
	energyScale_ = *(double*)(pbuffer); pbuffer += sizeof(double);
	rScale_ = *(double*)(pbuffer); pbuffer += sizeof(double);
	radialAngleScale_ = *(double*)(pbuffer); pbuffer += sizeof(double);
	azimutAngleScale_ = *(double*)(pbuffer); pbuffer += sizeof(double);
	noofParticles_ = *(unsigned*)(pbuffer); pbuffer += sizeof(unsigned);

	if (particles_ != nullptr)
		delete[] particles_;
	particles_ = new unsigned short[noofParticles_ * 4];
	memcpy(particles_, pbuffer, noofParticles_ * 4 * sizeof(unsigned short));

	// Настройка указателей
	for (unsigned i = 0; i < (unsigned)nThreads_; i++)
		// Начало самплинга в каждом потоке через равные промежутки
		currentParticleIdx_[i] = i * (noofParticles_ / nThreads_);
}

void mcSourceModelRadialDirect::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	//Источник задается в абсолютных координатах.
	//Поэтому ему не нужна матрица преобразования координат.

	int it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;

	os << "# Source: " << name_ << endl;

	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	// Концентрические круги
	double r = 0xffff * rScale_;
	for (it = 0; it < 360; it += da)
	{
		geomVector3D p = geomVector3D(r * sin(mPi * it), r * cos(mPi * it), z0_);
		os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
		p = geomVector3D(r * sin(mPi * (it + da)), r * cos(mPi * (it + da)), z0_);
		os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
		count++;
	}

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	for (it = 0; it < count; it++)
		os << "      " << 2 * it << ' ' << 2 * it + 1 << " -1" << endl;
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;
}

ostream& operator << (ostream& os, const mcSourceModelRadialDirect& s)
{
	os << (const mcSource&)s;
	os << "NAME = \t" << s.getName() << endl;
	os << "Z0 = \t" << s.z0_ << endl;
	os << "EnergyScale = \t" << s.energyScale_ << endl;
	os << "RadialAngleScale = \t" << s.radialAngleScale_ << endl;
	os << "AzimutAngleScale = \t" << s.azimutAngleScale_ << endl;
	os << "RScale = \t" << s.rScale_ << endl;
	os << "Noof particles = \t" << s.noofParticles_ << endl;
	os << endl;
	return os;
}
