#include "mcScoreSpectraFluenceSphere.h"
#include "mcThread.h"
#include "mcTransport.h"

mcScoreSpectraFluenceSphere::mcScoreSpectraFluenceSphere(const char* module_name, int nThreads, mc_particle_t pt, double ecut,
	int nr, int ne, double rmax, double emax, double radius)
	: mcScore(module_name, nThreads)
	, pt_(pt), nr_(nr), rstep_(rmax / (nr-1)), rmax_(rmax)
	, ne_(ne), estep_(emax / ne), emax_(emax), ecut_(ecut), radius_(radius)
	, energy_fluence_(nThreads), number_fluence_(nThreads), number_spectra_(nThreads), energy_spectra_(nThreads)
{
	setColor(0.5, 0.5, 0.0, 0.8);

	for (size_t i = 0; i < energy_fluence_.size(); i++)
	{
		energy_fluence_[i].resize(nr, 0);
		number_fluence_[i].resize(nr, 0);
		number_spectra_[i].resize(nr);
		energy_spectra_[i].resize(nr);
		for (int j = 0; j < nr; j++)
		{
			number_spectra_[i][j].resize(ne);
			energy_spectra_[i][j].resize(ne);
		}
	}
}

mcScoreSpectraFluenceSphere::~mcScoreSpectraFluenceSphere()
{
}

void mcScoreSpectraFluenceSphere::ScoreFluence(const mcParticle& particle)
{
	if (particle.t == pt_ && particle.ke >= ecut_)
	{
		int iThread = particle.thread_->id();
		double edep = particle.ke * particle.weight;
		etotal_[iThread] += edep;

		// !!! Scoring вызывается до перемещения частицы на поверхность объекта к которому привязана частица.
		// Нужно ее перенести здесь на поверхность.
		geomVector3D p = particle.p + (particle.u * (-particle.p.z() / particle.u.z()));

		// Положение частицы - это угол азимута градусах в предположении, 
		// что частица вылетает из центра координат.
		double angle = acos(particle.u.z()) * 180 / PI;
		int ridx = int(angle / rstep_ + 0.5);
		if (ridx >= 0 && ridx < nr_)
		{
			// Поток энергии
			energy_fluence_[iThread][ridx] += edep;
			number_fluence_[iThread][ridx] += particle.weight;

			int eidx = int(particle.ke / estep_);
			if (eidx < ne_)
			{
				// Спектр по количеству частиц
				number_spectra_[iThread][ridx][eidx] += particle.weight;

				// Спектр по энергии
				energy_spectra_[iThread][ridx][eidx] += edep;
			}
		}
	}
}

void mcScoreSpectraFluenceSphere::dumpVRML(ostream& os) const
{
	mcScore::dumpVRML(os);
	geomVector3D p = geomVector3D(0, 0, 0) * transport_->MT2W();

	os << "# Score: " << name_ << endl;
	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Sphere { radius " << radius_ << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}

void mcScoreSpectraFluenceSphere::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	int ir, ie;

	os << endl << "--------------------------------------------------------" << endl;
	os << "The integrated bremsstrahlung (SE0)" << endl;
	os << "--------------------------------------------------------" << endl;

	for (ir = 0; ir < nr_; ir++) os << "\t" << ir * rstep_;
	os << endl;
	os << "TE[MeV/sR]";
	for (ir = 0; ir < nr_; ir++)
		os << "\t" << (nFluenceBin(ir) / sAngle(ir));
	os << endl;

	os << endl << "--------------------------------------------------------" << endl;
	os << "BREMSSTRAHLUNG SPECTRA MEAN ENERGIES" << endl;
	os << "--------------------------------------------------------" << endl;

	for (ir = 0; ir < nr_; ir++) os << "\t" << ir * rstep_;
	os << endl;
	os << "ME[MeV]";
	for (ir = 0; ir < nr_; ir++)
		os << "\t" << (nFluenceBin(ir) > 0 ? eFluenceBin(ir) / nFluenceBin(ir) : 0.0);
	os << endl;

	os << endl << "---------------------------------------------------------------" << endl;
	os << "Azimatal number of particles spectra per sterradian per MeV [N]" << endl;
	os << "The differential bremsstrahlung yield (dS/dE)" << endl;
	os << "---------------------------------------------------------------" << endl;
	for (ir = 0; ir < nr_; ir++) os << "\t" << ir * rstep_;
	os << endl;
	for (ie = 0; ie < ne_; ie++)
	{
		os << (ie + 0.5) * estep_;
		for (ir = 0; ir < nr_; ir++)
			os << "\t" << (nSpectrumBin(ir, ie) / (sAngle(ir) * estep_));
		os << endl;
	}

	os << endl << "----------------------------------------------------" << endl;
	os << "Azimatal energy spectra per sterradian per MeV [MeV]" << endl;
	os << "(tally bremsstrablung spectra)" << endl;
	os << "----------------------------------------------------" << endl;
	for (ir = 0; ir < nr_; ir++) os << "\t" << ir * rstep_;
	os << endl;
	for (ie = 0; ie < ne_; ie++)
	{
		os << (ie + 0.5) * estep_;
		for (ir = 0; ir < nr_; ir++)
			os << "\t" << (eSpectrumBin(ir, ie) / (sAngle(ir) * estep_));
		os << endl;
	}
}

double mcScoreSpectraFluenceSphere::eFluenceBin(int ir) const
{
	double energy = 0;
	for (int i = 0; i < nThreads_; i++)
		energy += energy_fluence_[i][ir];
	return energy;
}

double mcScoreSpectraFluenceSphere::nFluenceBin(int ir) const
{
	double energy = 0;
	for (int i = 0; i < nThreads_; i++)
		energy += number_fluence_[i][ir];
	return energy;
}

double mcScoreSpectraFluenceSphere::nSpectrumBin(int ir, int ie) const
{
	double f = 0;
	for (int i = 0; i < nThreads_; i++)
		f += number_spectra_[i][ir][ie];
	return f;
}

double mcScoreSpectraFluenceSphere::eSpectrumBin(int ir, int ie) const
{
	double f = 0;
	for (int i = 0; i < nThreads_; i++)
		f += energy_spectra_[i][ir][ie];
	return f;
}

double mcScoreSpectraFluenceSphere::sAngle(int ir) const
{
	if (ir == 0)
		return 2 * PI * (1.0 - cos((0.5 * rstep_ * PI / 180.)));
	else
		return 2 * PI * (cos(((ir - 0.5) * rstep_) * PI / 180.) - cos(((ir + 0.5) * rstep_ * PI / 180.)));
}
