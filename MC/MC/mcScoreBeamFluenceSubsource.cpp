#include "mcScoreBeamFluenceSubsource.h"
#include "mcThread.h"

mcScoreBeamFluenceSubsource::mcScoreBeamFluenceSubsource(const char* name, int nThreads,
	enum mc_particle_t ptype, double zmin, double zmax, double focus,
	int nr, double rmax,
	int ne, double emax,
	int na, double amax)
	:ptype_(ptype), zmin_(zmin), zmax_(zmax), nThreads_(nThreads), focus_(focus), nr_(nr), rmax_(rmax)
	, ne_(ne), emax_(emax), na_(na), amax_(amax)
{
	if (name != nullptr)
		strcpy_s(name_, 32, name);

	rstep_ = rmax_ / nr_;
	estep_ = emax_ / ne_;
	astep_ = amax_ / na_;

	int len = nThreads_ * nr_;
	intencity_all_ = new double[len];
	memset(intencity_all_, 0, len * sizeof(double));
	intencity_ = new double*[nThreads];
	intencity_[0] = intencity_all_;
	for (int i = 1; i < nThreads_; i++)
		intencity_[i] = intencity_[i - 1] + nr_;

	len = nr_ * ne_;
	spectra_all_ = new double[nThreads_ * len];
	memset(spectra_all_, 0, nThreads_ * len * sizeof(double));
	spectra_ = new double*[nThreads];
	spectra_[0] = spectra_all_;
	for (int i = 1; i < nThreads_; i++)
		spectra_[i] = spectra_[i - 1] + len;

	len = nr_ * na_;
	angle_all_ = new double[nThreads_ * len];
	memset(angle_all_, 0, nThreads_ * len * sizeof(double));
	angle_ = new double*[nThreads];
	angle_[0] = angle_all_;
	for (int i = 1; i < nThreads_; i++)
		angle_[i] = angle_[i - 1] + len;
}

mcScoreBeamFluenceSubsource::~mcScoreBeamFluenceSubsource(void)
{
	if (intencity_all_) delete[] intencity_all_;
	if (intencity_) delete[] intencity_;

	if (spectra_all_) delete[] spectra_all_;
	if (spectra_) delete[] spectra_;

	if (angle_all_) delete[] angle_all_;
	if (angle_) delete[] angle_;
}

void mcScoreBeamFluenceSubsource::ScoreFluence(const mcParticle& particle, const geomVector3D& p)
{
	if (particle.t != ptype_ || particle.plast.z() < zmin_ || particle.plast.z() > zmax_)
		return;

	int iThread = particle.thread_->id();
	double edep = particle.ke * particle.weight;

	// Профиль потока энергии
	int ir = int(p.lengthXY() / rstep_);
	if (ir >= nr_)
		return;
	intencity_[iThread][ir] += edep;

	// Спектры в зависимости от радиуса
	int ie = int(particle.ke / estep_);
	if (ie < ne_)
		spectra_[iThread][ir*ne_ + ie] += edep;

	// Угловые разбросы в зависимости от радиуса
	geomVector3D vr(p.x(), p.y(), focus_);
	vr.normalize();
	double a = acos(vr * particle.u) * 180 / PI;
	int ia = int(a / astep_);
	if (ia < na_)
		angle_[iThread][ir*na_ + ia] += edep;
}

double mcScoreBeamFluenceSubsource::Intencity(int iThread, int ir) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluenceSubsource::Intencity: thread exceed container size");
	else if (ir < 0 || ir >= nr_)
		return 0;
	return intencity_[iThread][ir] / (PI * rstep_ * rstep_ * (2 * ir + 1));
}

double mcScoreBeamFluenceSubsource::Intencity(int ir) const
{
	double f = 0;
	if (ir >= 0 && ir < nr_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += intencity_[i][ir];
		f /= PI * rstep_ * rstep_ * (2 * ir + 1);
	}
	return f;
}

double mcScoreBeamFluenceSubsource::Spectrum(int iThread, int ir, int ie) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluenceSubsource::Spectrum: thread exceed container size");
	else if (ir < 0 || ir >= nr_ || ie < 0 || ie >= ne_)
		return 0;
	return spectra_[iThread][ir*ne_ + ie];
}

double mcScoreBeamFluenceSubsource::Spectrum(int ir, int ie) const
{
	double f = 0;
	if (ir >= 0 && ir < nr_)
	{
		for (int i = 0; i < nThreads_; i++)
			//f += spectra_[i][ir*ne_ + ie];
			f += Spectrum(i, ir, ie);
	}
	return f;
}

double mcScoreBeamFluenceSubsource::SpectrumSigma(int ir, int ie) const
{
	double f1 = 0, f2 = 0;
	if (ir >= 0 && ir < nr_)
	{
		for (int i = 0; i < nThreads_; i++)
		{
			double f = Spectrum(i, ir, ie);
			f1 += f;
			f2 += f * f;
		}
	}
	return sqrt((f2 - f1*f1 / nThreads_) / nThreads_);
}

double mcScoreBeamFluenceSubsource::Angle(int iThread, int ir, int ia) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluenceSubsource::Spectrum: thread exceed container size");
	else if (ir < 0 || ir >= nr_ || ia < 0 || ia >= na_)
		return 0;
	return angle_[iThread][ir*na_ + ia];
}

double mcScoreBeamFluenceSubsource::Angle(int ir, int ia) const
{
	double f = 0;
	if (ir >= 0 && ir < nr_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += Angle(i, ir, ia);
	}
	return f;
}

double mcScoreBeamFluenceSubsource::AngleSigma(int ir, int ia) const
{
	double f1 = 0, f2 = 0;
	if (ir >= 0 && ir < nr_)
	{
		for (int i = 0; i < nThreads_; i++)
		{
			double f = Angle(i, ir, ia);
			f1 += f;
			f2 += f * f;
		}
	}
	return sqrt((f2 - f1*f1 / nThreads_) / nThreads_);
}

void mcScoreBeamFluenceSubsource::dumpStatistic(ostream& os) const
{
	int ir, ie, ia;

	os << endl;
	os << "Subsource name: " << name_ << endl;
	os << "---------------------------------------" << endl;

	os << endl;
	os << "Fluence intencity profile" << endl;
	os << endl;

	for (ir = 0; ir < nr_; ir++)
		os << '\t' << (double(ir) + 0.5) * rstep_;
	os << endl;

	for (ir = 0; ir < nr_; ir++)
		os << '\t' << Intencity(ir);
	os << endl;

	os << endl;
	os << "Spectra versus radius" << endl;
	os << endl;

	double* spec_tot = new double[ne_];
	double* sigma_tot = new double[ne_];
	for (ie = 0; ie < ne_; ie++)
	{
		spec_tot[ie] = 0;
		sigma_tot[ie] = 0;
		os << '\t' << (double(ie) + 0.5) * estep_;
	}
	for (ie = 0; ie < ne_; ie++)
		os << '\t' << (double(ie) + 0.5) * estep_;
	os << endl;

	double ett = 0;
	for (ir = 0; ir < nr_; ir++)
	{
		// Нормируем спектры таким образом, что суммарное кол-во частиц равно 1.
		os << (double(ir) + 0.5) * rstep_;

		double et = 0;
		for (ie = 0; ie < ne_; ie++)
		{
			double f = Spectrum(ir, ie);
			spec_tot[ie] += f;
			et += f;
			ett += f;
			f = SpectrumSigma(ir, ie);
			sigma_tot[ie] += f * f;
		}

		for (ie = 0; ie < ne_; ie++)
		{
			if (et == 0)
				os << "\t0";
			else
				os << '\t' << Spectrum(ir, ie) / et;
		}

		// Sigma
		for (ie = 0; ie < ne_; ie++)
		{
			if (et == 0)
				os << "\t0";
			else
				os << '\t' << SpectrumSigma(ir, ie) / et;
		}
		os << endl;
	}

	// Суммарный спектр
	os << endl;
	os << "Суммарный спектр" << endl;
	os << endl;

	for (ie = 0; ie < ne_; ie++)
		os << '\t' << (double(ie) + 0.5) * estep_;
	for (ie = 0; ie < ne_; ie++)
		os << '\t' << (double(ie) + 0.5) * estep_;
	os << endl;

	for (ie = 0; ie < ne_; ie++)
		os << '\t' << spec_tot[ie] / ett;
	for (ie = 0; ie < ne_; ie++)
		os << '\t' << sqrt(sigma_tot[ie]) / ett;

	os << endl;
	delete[] spec_tot;
	delete[] sigma_tot;

	os << endl;
	os << "Angle spread versus radius (per solid angle)" << endl;
	os << endl;

	double* ang_tot = new double[na_];
	double* ang_sigma_tot = new double[na_];
	for (ia = 0; ia < na_; ia++)
	{
		ang_tot[ia] = 0;
		ang_sigma_tot[ia] = 0;
		os << '\t' << (double(ia) + 0.5) * astep_;
	}
	for (ia = 0; ia < na_; ia++)
		os << '\t' << (double(ia) + 0.5) * astep_;
	os << endl;

	double att = 0;
	for (ir = 0; ir < nr_; ir++)
	{
		os << (double(ir) + 0.5) * rstep_;

		double at = 0;
		for (ia = 0; ia < na_; ia++)
		{
			double f = Angle(ir, ia) / (cos(double(ia + 1)*PI / 180) - cos(double(ia)*PI / 180));
			ang_tot[ia] += f;
			at += f;
			att += f;
			f = AngleSigma(ir, ia) / (cos(double(ia + 1)*PI / 180) - cos(double(ia)*PI / 180));
			ang_sigma_tot[ia] += f * f;
		}

		for (ia = 0; ia < na_; ia++)
		{
			if (at == 0)
				os << "\t0";
			else
				os << '\t' << Angle(ir, ia) / ((cos(double(ia + 1)*PI / 180) - cos(double(ia)*PI / 180)) * at);
		}

		// Sigma
		for (ia = 0; ia < na_; ia++)
		{
			if (at == 0)
				os << "\t0";
			else
				os << '\t' << AngleSigma(ir, ia) / ((cos(double(ia + 1)*PI / 180) - cos(double(ia)*PI / 180)) * at);
		}
		os << endl;
	}

	// Суммарное угловое распределение
	os << endl;
	os << "Суммарное угловое распределение" << endl;
	os << endl;

	for (ia = 0; ia < na_; ia++)
		os << '\t' << (double(ia) + 0.5) * astep_;
	for (ia = 0; ia < na_; ia++)
		os << '\t' << (double(ia) + 0.5) * astep_;
	os << endl;

	for (ia = 0; ia < na_; ia++)
		os << '\t' << ang_tot[ia] / att;
	for (ia = 0; ia < na_; ia++)
		os << '\t' << sqrt(ang_sigma_tot[ia]) / att;
	os << endl;
	delete[] ang_tot;
	delete[] ang_sigma_tot;
}
