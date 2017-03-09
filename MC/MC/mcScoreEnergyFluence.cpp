#include "mcScoreEnergyFluence.h"
#include "mcThread.h"
#include "mcTransport.h"

mcScoreEnergyFluence::mcScoreEnergyFluence(const char* module_name, int nThreads, mc_particle_t pt, int nr, double rmax)
	:mcScore(module_name, nThreads), pt_(pt), nr_(nr), rstep_(rmax / nr), rmax_(rmax), fluence_(nThreads)
{
	for (size_t i = 0; i < fluence_.size(); i++)
		fluence_[i].resize(nr, 0);
}

mcScoreEnergyFluence::~mcScoreEnergyFluence()
{
}

void mcScoreEnergyFluence::ScoreFluence(const mcParticle& particle)
{
	if (particle.t == pt_)
	{
		int iThread = particle.thread_->id();
		double edep = particle.ke * particle.weight;
		etotal_[iThread] += edep;

		// !!! Scoring вызывается до перемещения частицы на поверхность объекта к которому привязана частица.
		// Нужно ее перенести здесь на поверхность.
		geomVector3D p = particle.p + (particle.u * (-particle.p.z() / particle.u.z()));

		double r = p.lengthXY();
		int ridx = int(r / rstep_);
		if (ridx < nr_)
			fluence_[iThread][ridx] += edep;
	}
}

void mcScoreEnergyFluence::dumpVRML(ostream& os) const
{
	os << "# mcScoreEnergyFluence Score: " << name_ << endl;
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

void mcScoreEnergyFluence::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	double sf = 1 / (2 * PI * rstep_ * rstep_);

	os << endl << "Радиально симметричные профили потока энергии" << endl;
	os << "----------------------------------------------------- " << endl;
	int i, ni = (int)fluence_[0].size();

	for (i = -ni + 1; i <= 0; i++)
	{
		double energy = 0;
		for (unsigned j = 0; j < fluence_.size(); j++)
			energy += fluence_[j][-i];

		os << (i - 0.5) * rstep_ << "\t" << (energy * sf / (-2 * i + 1)) << endl;
	}

	for (i = 0; i < ni; i++)
	{
		double energy = 0;
		for (unsigned j = 0; j < fluence_.size(); j++)
			energy += fluence_[j][i];

		os << (i + 0.5) * rstep_ << "\t" << (energy * sf / (2 * i + 1)) << endl;
	}
}
