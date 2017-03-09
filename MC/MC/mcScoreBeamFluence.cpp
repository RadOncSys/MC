#include "mcScoreBeamFluence.h"
#include "mcParticle.h"
#include "mcThread.h"
#include "mcTransport.h"

mcScoreBeamFluence::mcScoreBeamFluence(const char* module_name, int nThreads, int nr, double rmax)
	:mcScore(module_name, nThreads)
	, nsub_(0), nr_(nr), rmax_(rmax)
{
	rstep_ = rmax_ / nr_;
}

mcScoreBeamFluence::~mcScoreBeamFluence()
{
	for (int i = 0; i < nsub_; i++)
		delete subsources_[i];
}

void mcScoreBeamFluence::ScoreFluence(const mcParticle& particle)
{
	int iThread = particle.thread_->id();
	double edep = particle.ke * particle.weight;
	etotal_[iThread] += edep;

	// !!! Scoring вызываетс€ до перемещени€ частицы на поверхность объекта к которому прив€зана частица.
	// Ќужно ее перенести здесь на поверхность.
	geomVector3D p = particle.p + (particle.u * (-particle.p.z() / particle.u.z()));

	for (int i = 0; i < nsub_; i++)
		subsources_[i]->ScoreFluence(particle, p);
}

void mcScoreBeamFluence::addSubsource(mcScoreBeamFluenceSubsource* s)
{
	if (nsub_ == 10)
		throw std::exception("mcScoreBeamFluence::addSubsource: too many subsources");
	subsources_[nsub_++] = s;
}

void mcScoreBeamFluence::dumpVRML(ostream& os) const
{
	os << "# mcScoreBeamFluence Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int ir, it, count = 0;
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

	//  онцентрические круги
	for (ir = 1; ir <= nr_; ir++) {
		double r = rstep_ * ir;
		for (it = 0; it < 360; it += da) {
			geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), 0) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r*sin(mPi*(it + da)), r*cos(mPi*(it + da)), 0) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			count++;
		}
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

void mcScoreBeamFluence::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);
	for (int i = 0; i < nsub_; i++)
		subsources_[i]->dumpStatistic(os);
}
