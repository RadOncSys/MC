#include "mcScorePhaseSpaceDirect.h"
#include "mcParticle.h"
#include "mcThread.h"
#include "mcTransport.h"
#include "mcSourceModelRadialDirect.h"

#include <vector>
using namespace std;

#define MAXNOOF_PARTICLES_PER_THREAD	10000000

mcScorePhaseSpaceDirect::mcScorePhaseSpaceDirect(
	const char* module_name, int nThreads, enum mc_particle_t ptype,
	double emax, double rmax, const char* stat_file)
	:mcScore(module_name, nThreads)
	, ptype_(ptype)
	, emax_(emax)
	, rmax_(rmax)
	, model_file_(stat_file)
{
	data_ = new unsigned short [MAXNOOF_PARTICLES_PER_THREAD * nThreads * 4];
	memset(data_, 0, MAXNOOF_PARTICLES_PER_THREAD * nThreads * 4 * sizeof(unsigned short));
	threadIndexes_.resize(nThreads, 0);

	energyScale_ = emax_ / 0xffff;
	rScale_ = rmax_ / 0xffff;
	radialAngleScale_ = 2.0 / 0xffff;
	azimutAngleScale_ = 2.0 / 0xffff;

	// Debug
	nr_ = 50;
	dr_ = rmax_ * 1.2 / nr_;
	weight_.resize(nr_, 0);
	fluence_.resize(nr_, 0);
	rangle_.resize(nr_, 0);
	aangle_.resize(nr_, 0);
	rangle2_.resize(nr_, 0);
	aangle2_.resize(nr_, 0);
}

mcScorePhaseSpaceDirect::~mcScorePhaseSpaceDirect()
{
	if (data_ != nullptr)
		delete [] data_;
}

void mcScorePhaseSpaceDirect::ScoreFluence(const mcParticle& particle)
{
	if (particle.t != ptype_ || particle.u.z() <= 0)
		return;

	int iThread = particle.thread_->id();
	etotal_[iThread] += particle.ke * particle.weight;

	unsigned pidx = threadIndexes_[iThread];

	if (pidx < MAXNOOF_PARTICLES_PER_THREAD)
	{
		unsigned short* sp = data_ + (iThread * MAXNOOF_PARTICLES_PER_THREAD + pidx) * 4;

		double x = particle.p.x(), y = particle.p.y();
		double r = sqrt(x * x + y * y);
		if (r < rmax_)
		{
			sp[0] = (unsigned short)(particle.ke / energyScale_);
			sp[1] = (unsigned short)(r / rScale_);

			double sinPhi = y / r, cosPhi = x / r;

			double ux = particle.u.x(), uy = particle.u.y();
			double xa = ux * cosPhi + uy * sinPhi, ya = -ux * sinPhi + uy * cosPhi;
			sp[2] = (unsigned short)((xa + 1.0) / radialAngleScale_);
			sp[3] = (unsigned short)((ya + 1.0) / azimutAngleScale_);

			threadIndexes_[iThread]++;

			// Debug
			unsigned ir = unsigned(r / dr_);
			weight_[ir] += particle.weight;
			fluence_[ir] += particle.ke * particle.weight;
			rangle_[ir] += xa;
			rangle2_[ir] += xa * xa;
			aangle_[ir] += ya;
			aangle2_[ir] += ya * ya;
		}
	}
}

void mcScorePhaseSpaceDirect::dumpStatistic(ostream& os) const
{
	// Создаем файл модели
	auto source = std::make_unique<mcSourceModelRadialDirect>("Radial symmetric particle colleection source", 1, 0);
	source->Init(emax_, rmax_);
	source->SetParticles(data_, MAXNOOF_PARTICLES_PER_THREAD, threadIndexes_);

	int len;
	void* memobj = source->saveToMemory(len);
	if (memobj == nullptr)
		throw std::exception("mcScorePhaseSpaceDirect:: Cannot create model object in memory");

	FILE* modelfile;
	if (fopen_s(&modelfile, model_file_.c_str(), "wb") != 0)
		throw std::exception("mcScorePhaseSpaceDirect:: Cannot open model file for writing");
	fwrite(memobj, 1, len, modelfile);
	fclose(modelfile);

	free(memobj);

	mcScore::dumpStatistic(os);

	os << "ptype_:\t" << ptype_ << endl;
	os << "emax_:\t" << emax_ << endl;
	os << "rmax_:\t" << rmax_ << endl;
	os << "model_file_:\t" << model_file_.c_str() << endl;

	os << endl << "Source model";
	os << endl << "============" << endl;
	os << *source;

	// Интегральные характеристики для перехвата грубых ошибок

	auto np = source->GetNoofParticles();
	const unsigned short* particles = source->GetParticles();

	unsigned nr = 50, badcount = 0;
	double dr = rmax_ * 1.2 / nr;
	std::vector<double> weight(nr, 0);
	std::vector<double> fluence(nr, 0);
	std::vector<double> rangle(nr, 0);
	std::vector<double> aangle(nr, 0);
	std::vector<double> rangle2(nr, 0);
	std::vector<double> aangle2(nr, 0);

	for (unsigned i = 0; i < np; i++)
	{
		const unsigned short* sp = particles + i * 4;
		double r = sp[1] * rScale_;
		unsigned ir = unsigned(r / dr);
		if (ir < nr)
		{
			double xa = sp[2] * radialAngleScale_ - 1.0, ya = sp[3] * azimutAngleScale_ - 1.0;
			weight[ir] += 1.0;
			fluence[ir] += sp[0] * energyScale_;
			rangle[ir] += xa;
			rangle2[ir] += xa * xa;
			aangle[ir] += ya;
			aangle2[ir] += ya * ya;
		}
		else badcount++;
	}

	// 3-й набор - это уже самплинг из модели
	unsigned nsplit = 10;
	unsigned nsamp = nsplit * source->GetNoofParticles();
	source->SetSplitting(nsplit);
	mcParticle particle;
	mcThread thread;
	thread.setId(0);
	std::vector<double> s_weight(nr, 0);
	std::vector<double> s_fluence(nr, 0);
	std::vector<double> s_rangle(nr, 0);
	std::vector<double> s_aangle(nr, 0);
	std::vector<double> s_rangle2(nr, 0);
	std::vector<double> s_aangle2(nr, 0);

	for (unsigned i = 0; i < nsamp; i++)
	{
		source->sample(particle, &thread);

		double x = particle.p.x(), y = particle.p.y();
		double r = sqrt(x * x + y * y);
		if (r < rmax_)
		{
			double sinPhi = y / r, cosPhi = x / r;
			double ux = particle.u.x(), uy = particle.u.y();
			double xa = ux * cosPhi + uy * sinPhi, ya = -ux * sinPhi + uy * cosPhi;

			unsigned ir = unsigned(r / dr_);
			s_weight[ir] += particle.weight;
			s_fluence[ir] += particle.ke * particle.weight;
			s_rangle[ir] += xa * particle.weight;
			s_rangle2[ir] += xa * xa * particle.weight;
			s_aangle[ir] += ya * particle.weight;
			s_aangle2[ir] += ya * ya * particle.weight;
		}
	}

	os << endl << endl << "Bad count = " << badcount << endl << endl;

	os << endl << "Energy fluence as function of radius";
	os << endl << "====================================" << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (dr * ( i + 0.5));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (fluence[i] / (PI * dr * dr * (2 * i + 1)));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (fluence_[i] / (PI * dr_ * dr_ * (2 * i + 1)));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (s_fluence[i] / (PI * dr * dr * (2 * i + 1)));
	os << endl;

	os << endl << "Mean energy as function of radius";
	os << endl << "=================================" << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (dr * (i + 0.5));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight[i] > 0 ? fluence[i] / weight[i] : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight_[i] > 0 ? fluence_[i] / weight_[i] : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (s_weight[i] > 0 ? s_fluence[i] / s_weight[i] : 0);
	os << endl;

	os << endl << "Mean radial angle";
	os << endl << "=================" << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (dr * (i + 0.5));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight[i] > 0 ? rangle[i] / weight[i] : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight_[i] > 0 ? rangle_[i] / weight_[i] : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (s_weight[i] > 0 ? s_rangle[i] / s_weight[i] : 0);
	os << endl;

	os << endl << "Mean azimutal angle";
	os << endl << "===================" << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (dr * (i + 0.5));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight[i] > 0 ? aangle[i] / weight[i] : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight_[i] > 0 ? aangle_[i] / weight_[i] : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (s_weight[i] > 0 ? s_aangle[i] / s_weight[i] : 0);
	os << endl;

	os << endl << "STDEV radial angle";
	os << endl << "==================" << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (dr * (i + 0.5));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight[i] > 0 ? sqrt((rangle2[i] - (rangle[i] * rangle[i]) / (weight[i] * weight[i])) / weight[i]) : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight_[i] > 0 ? sqrt((rangle2_[i] - (rangle_[i] * rangle_[i]) / (weight_[i] * weight_[i])) / weight_[i]) : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (s_weight[i] > 0 ? sqrt((s_rangle2[i] - (s_rangle[i] * s_rangle[i]) / (s_weight[i] * s_weight[i])) / s_weight[i]) : 0);
	os << endl;

	os << endl << "STDEV azimutal angle";
	os << endl << "====================" << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (dr * (i + 0.5));
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight[i] > 0 ? sqrt((aangle2[i] - (aangle[i] * aangle[i]) / (weight[i] * weight[i])) / weight[i]) : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (weight_[i] > 0 ? sqrt((aangle2_[i] - (rangle_[i] * rangle_[i]) / (weight_[i] * weight_[i])) / weight_[i]) : 0);
	os << endl;
	for (unsigned i = 0; i < nr; i++) os << "\t" << (s_weight[i] > 0 ? sqrt((s_aangle2[i] - (s_aangle[i] * s_aangle[i]) / (s_weight[i] * s_weight[i])) / s_weight[i]) : 0);
	os << endl;
}

void mcScorePhaseSpaceDirect::dumpVRML(ostream& os) const
{
	os << "# mcScorePhaseSpaceDirect Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;

	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	// Круг
	double r = rmax_;
	for (it = 0; it < 360; it += da) {
		geomVector3D p = geomVector3D(r * sin(mPi * it), r * cos(mPi * it), 0) * mttow;
		os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
		p = geomVector3D(r * sin(mPi * (it + da)), r * cos(mPi * (it + da)), 0) * mttow;
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
