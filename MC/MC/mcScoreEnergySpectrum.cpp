#include "mcScoreEnergySpectrum.h"
#include "mcThread.h"
#include "mcTransport.h"

mcScoreEnergySpectrum::mcScoreEnergySpectrum(const char* module_name, int nThreads, mc_particle_t pt, int ne, double emax, double rmax)
	:mcScore(module_name, nThreads), pt_(pt), ne_(ne), emax_(emax), estep_(emax / ne), rmax_(rmax * rmax), espec_(nThreads)
{
	for (size_t i = 0; i < espec_.size(); i++)
		espec_[i].resize(ne, 0);
}

mcScoreEnergySpectrum::~mcScoreEnergySpectrum()
{
}

void mcScoreEnergySpectrum::ScoreFluence(const mcParticle& particle)
{
	if (particle.t == pt_ && particle.p.sqLengthXY() < rmax_)
	{
		int iThread = particle.thread_->id();
		double edep = particle.ke * particle.weight;
		etotal_[iThread] += edep;

		int eidx = int(particle.ke / estep_);
		if (eidx < ne_)
			espec_[iThread][eidx] += particle.weight;
	}
}

void mcScoreEnergySpectrum::dumpVRML(ostream& os) const
{
	os << "# mcScoreEnergySpectrum Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;
	double r = sqrt(rmax_);

	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	// Концентрический круг
	for (it = 0; it < 360; it += da) {
		geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), 0) * mttow;
		os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
		p = geomVector3D(r*sin(mPi*(it + da)), r*cos(mPi*(it + da)), 0) * mttow;
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

void mcScoreEnergySpectrum::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	os << endl << "Распределение по энергии (МэВ)" << endl;
	os << "--------------------------------------" << endl;
	for (unsigned i = 0; i < espec_[0].size(); i++)
	{
		double energy = 0;
		for (unsigned j = 0; j < espec_.size(); j++)
			energy += espec_[j][i];

		os << (i + 0.5) * estep_ << "\t" << energy << endl;
	}
}
