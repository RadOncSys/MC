#include "mcScoreSpectraFluence.h"
#include "mcThread.h"
#include "mcTransport.h"

mcScoreSpectraFluence::mcScoreSpectraFluence(const char* module_name, int nThreads, mc_particle_t pt, double ecut, 
	int nr, int ne, double rmax, double emax)
	: mcScore(module_name, nThreads)
	, pt_(pt), nr_(nr), rstep_(rmax / nr), rmax_(rmax)
	, ne_(ne), estep_(emax / ne), emax_(emax), ecut_(ecut)
	, energy_fluence_(nThreads), number_spectra_(nThreads), energy_spectra_(nThreads)
{
	for (size_t i = 0; i < energy_fluence_.size(); i++)
	{
		energy_fluence_[i].resize(nr, 0);
		number_spectra_[i].resize(nr);
		energy_spectra_[i].resize(nr);
		for (int j = 0; j < nr; j++)
		{
			number_spectra_[i][j].resize(ne);
			energy_spectra_[i][j].resize(ne);
		}
	}
}

mcScoreSpectraFluence::~mcScoreSpectraFluence()
{
}

void mcScoreSpectraFluence::ScoreFluence(const mcParticle& particle)
{
	if (particle.t == pt_ && particle.ke >= ecut_)
	{
		int iThread = particle.thread_->id();
		double edep = particle.ke * particle.weight;
		etotal_[iThread] += edep;

		// !!! Scoring вызывается до перемещения частицы на поверхность объекта к которому привязана частица.
		// Нужно ее перенести здесь на поверхность.
		geomVector3D p = particle.p + (particle.u * (-particle.p.z() / particle.u.z()));

		double r = p.lengthXY();
		int ridx = int(r / rstep_);
		if (ridx >= nr_ || ridx < 0) return;

		// Поток энергии
		energy_fluence_[iThread][ridx] += edep;

		// Спектр по количеству частиц
		int eidx = int(particle.ke / estep_);
		if (eidx >= ne_) return;

		number_spectra_[iThread][ridx][eidx] += particle.weight;

		// Спектр по энергии
		energy_spectra_[iThread][ridx][eidx] += edep;
	}
}

void mcScoreSpectraFluence::dumpVRML(ostream& os) const
{
	os << "# mcScoreSpectraFluence Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;
	double r = rmax_;

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
	for (int ir = 1; ir <= nr_; ir++) {
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

void mcScoreSpectraFluence::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	int ir, ie;
	double sf = 1 / (2 * PI * rstep_ * rstep_);

	os << endl << "Радиально симметричный профиль потока энергии";
	os << endl << "---------------------------------------------" << endl;

	for (ir = -nr_ + 1; ir <= 0; ir++) os << (ir - 0.5) * rstep_ << "\t";
	for (ir = 0; ir < nr_; ir++) os << (ir + 0.5) * rstep_ << "\t";
	os << endl;

	for (ir = -nr_ + 1; ir <= 0; ir++) os << (eFluenceBin(-ir) * sf / (-2 * ir + 1)) << "\t";
	for (ir = 0; ir < nr_; ir++) os << (eFluenceBin(ir) * sf / (2 * ir + 1)) << "\t";
	os << endl;

	os << endl << "Радиальные спектры по количеству частиц";
	os << endl << "---------------------------------------" << endl;
	for (ie = 0; ie < ne_; ie++) os << "\t" << (ie - 0.5) * estep_;
	os << endl;
	for (ir = 0; ir < nr_; ir++)
	{
		os << (ir + 0.5) * rstep_;
		for (ie = 0; ie < ne_; ie++) os << "\t" << (nSpectrumBin(ir, ie) * sf / (2 * ir + 1));
		os << endl;
	}

	os << endl << "Радиальные спектры по энергии частиц";
	os << endl << "------------------------------------" << endl;
	for (ie = 0; ie < ne_; ie++) os << "\t" << (ie - 0.5) * estep_;
	os << endl;
	for (ir = 0; ir < nr_; ir++)
	{
		os << (ir + 0.5) * rstep_;
		for (ie = 0; ie < ne_; ie++) os << "\t" << (eSpectrumBin(ir, ie) * sf / (2 * ir + 1));
		os << endl;
	}
}

double mcScoreSpectraFluence::eFluenceBin(int ir) const
{
	double energy = 0;
	for (int i = 0; i < nThreads_; i++)
		energy += energy_fluence_[i][ir];
	return energy;
}

double mcScoreSpectraFluence::nSpectrumBin(int ir, int ie) const
{
	double f = 0;
	for (int i = 0; i < nThreads_; i++)
		f += number_spectra_[i][ir][ie];
	return f;
}

double mcScoreSpectraFluence::eSpectrumBin(int ir, int ie) const
{
	double f = 0;
	for (int i = 0; i < nThreads_; i++)
		f += energy_spectra_[i][ir][ie];
	return f;
}
