#include "mcScoreAcceleratedBeam.h"
#include "mcParticle.h"

mcScoreAcceleratedBeam::mcScoreAcceleratedBeam(const char* module_name, int nThreads, enum mc_particle_t ptype,
	int ne, int nr, int naxial, int nthet, double emin, double emax, double rmax, double thetmax)
	:mcScore(module_name, nThreads)
	, ptype_(ptype)
	, ne_(ne), nr_(nr), nthet_(nthet)
	, naxial_(naxial), emin_(emin)
	, emax_(emax)
	, rmax_(rmax)
	, skipped_energy_(0)
{
	if (nThreads != 1)
		throw std::exception("mcScoreAcceleratedBeam: данный скоринг поддерживает только один поток");

	de_ = (emax_ - emin_) / ne_;
	dr_ = rmax_ / nr_;
	thetmax_ = thetmax * PI / 180.0;
	dthet_ = thetmax_ / nthet_;
	daxial_ = 2 * PI / naxial_;

	int i;
	espec_.resize(ne_, 0);
	rspec_.resize(nr_, 0);
	aspec_.resize(naxial_, 0);
	tspec_.resize(nthet_, 0);

	evr_.resize(nr_);
	for (i = 0; i < (int)evr_.size(); i++) evr_[i].resize(ne_, 0);

	eva_.resize(naxial_);
	for (i = 0; i < (int)eva_.size(); i++) eva_[i].resize(ne_, 0);

	tva_.resize(naxial_);
	for (i = 0; i < (int)tva_.size(); i++) tva_[i].resize(nthet_, 0);

	wvxy_.resize(nr_ * 2);
	for (i = 0; i < (int)wvxy_.size(); i++) wvxy_[i].resize(nr_ * 2, 0);

	wvvxvy_.resize(nthet_ * 2);
	for (i = 0; i < (int)wvvxvy_.size(); i++) wvvxvy_[i].resize(nthet_ * 2, 0);
}

void mcScoreAcceleratedBeam::ScoreFluence(const mcParticle& particle)
{
	if (particle.t != ptype_)
		return;

	// Координаты
	double e = particle.ke;
	double w = e * particle.weight;
	double x = particle.p(0), y = particle.p(1);
	double r = sqrt(x*x + y*y);
	double vx = particle.u(0), vy = particle.u(1), vz = particle.u(2);

	etotal_[0] += w;

	// Рассчитываем индексы
	int ie = int((particle.ke - emin_) / de_);
	if (ie < 0 || ie >= ne_) { skipped_energy_ += w; return; }

	int ir = int(r / dr_);
	int ix = int(nr_ + x / dr_);
	int iy = int(nr_ + y / dr_);
	if (ir >= nr_) { skipped_energy_ += w; return; };
	if (ix < 0 || ix >= 2 * nr_) { skipped_energy_ += w; return; };
	if (iy < 0 || iy >= 2 * nr_) { skipped_energy_ += w; return; };

	double fa = x > 0 ? atan(y / x) : x < 0 ? (PI + atan(y / x)) : y > 0 ? 0.5 * PI : 1.5 * PI;
	if (fa < 0) fa += 2 * PI;
	int ia = int(fa / daxial_);
	if (ia < 0 || ia >= naxial_) { skipped_energy_ += w; return; };

	double fva = vx > 0 ? atan(vy / vx) : vx < 0 ? (PI + atan(vy / vx)) : vy > 0 ? 0.5 * PI : 1.5 * PI;
	if (fva < 0) fva += 2 * PI;
	int iva = int(fva / daxial_);
	if (iva < 0 || iva >= naxial_) { skipped_energy_ += w; return; };

	int itheta = int(acos(vz) / dthet_);
	if (itheta >= nthet_) { skipped_energy_ += w; return; };

	int ivx = int(nthet_ + asin(vx) / dthet_);
	int ivy = int(nthet_ + asin(vy) / dthet_);
	if (ivx < 0 || ivx >= 2 * nthet_) { skipped_energy_ += w; return; };
	if (ivy < 0 || ivy >= 2 * nthet_) { skipped_energy_ += w; return; };

	// Наполняем массивы
	espec_[ie] += w;
	rspec_[ir] += w;
	aspec_[iva] += w;
	tspec_[itheta] += w;

	evr_[ir][ie] += w;
	eva_[iva][ie] += w;
	tva_[iva][itheta] += w;
	wvxy_[iy][ix] += w;
	wvvxvy_[ivy][ivx] += w;
}

void mcScoreAcceleratedBeam::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	os << "Распределение частиц, падающих на радиационную мишень электронного ускорителя" << endl << endl;

	os << "ptype_:\t" << ptype_ << endl;
	os << "ne_:\t" << ne_ << endl;
	os << "nr_:\t" << nr_ << endl;
	os << "nthet_:\t" << nthet_ << endl;
	os << "naxial_:\t" << naxial_ << endl;
	os << "emin_:\t" << emin_ << endl;
	os << "emax_:\t" << emax_ << endl;
	os << "rmax_:\t" << rmax_ << endl;
	os << "thetmax_:\t" << thetmax_ << endl;
	os << endl;
	os << "Суммарная энергия:\t" << etotal_[0] << endl;
	os << "Пропущенная энергия:\t" << skipped_energy_ << endl;

	int i, j;

	os << endl << "Суммарный энергетический спектр" << endl << endl;
	for (i = 0; i < (int)espec_.size(); i++)
		os << "\t" << emin_ + (i + 0.5) * de_;
	os << endl;
	for (i = 0; i < (int)espec_.size(); i++)
		os << "\t" << espec_[i];
	os << endl;

	os << endl << "Распределение энергии по радиусу" << endl << endl;
	for (i = 0; i < (int)rspec_.size(); i++)
		os << "\t" << (i + 0.5) * dr_;
	os << endl;
	for (i = 0; i < (int)rspec_.size(); i++)
		os << "\t" << rspec_[i] / (2 * i + 1);
	os << endl;

	os << endl << "Распределение энергии по азимуту положения" << endl << endl;
	for (i = 0; i < (int)aspec_.size(); i++)
		os << "\t" << (i + 0.5) * daxial_ * 180 / PI;
	os << endl;
	for (i = 0; i < (int)aspec_.size(); i++)
		os << "\t" << aspec_[i];
	os << endl;

	os << endl << "Распределение энергии по направлению" << endl << endl;
	for (i = 0; i < (int)tspec_.size(); i++)
		os << "\t" << (i + 0.5) * dthet_ * 180 / PI;
	os << endl;
	for (i = 0; i < (int)tspec_.size(); i++)
		os << "\t" << tspec_[i] / (2 * i + 1);
	os << endl;

	os << endl << "Энергетические спектры в зависимости от радиуса" << endl << endl;
	for (i = 0; i < (int)espec_.size(); i++)
		os << "\t" << emin_ + (i + 0.5) * de_;
	os << endl;
	for (i = 0; i < (int)evr_.size(); i++)
	{
		os << (i + 0.5) * dr_;
		for (j = 0; j < (int)evr_[0].size(); j++)
			os << "\t" << evr_[i][j];
		os << endl;
	}

	os << endl << "Энергетические спектры в зависимости от азимута" << endl << endl;
	for (i = 0; i < (int)aspec_.size(); i++)
		os << "\t" << emin_ + (i + 0.5) * de_;
	os << endl;
	for (i = 0; i < (int)eva_.size(); i++)
	{
		os << (i + 0.5) * daxial_ * 180 / PI;
		for (j = 0; j < (int)eva_[0].size(); j++)
			os << "\t" << eva_[i][j];
		os << endl;
	}

	os << endl << "Угловые распределения в зависимости от азимута" << endl << endl;
	for (i = 0; i < (int)tspec_.size(); i++)
		os << "\t" << (i + 0.5) * dthet_ * 180 / PI;
	os << endl;
	for (i = 0; i < (int)tva_.size(); i++)
	{
		os << (i + 0.5) * daxial_ * 180 / PI;
		for (j = 0; j < (int)tva_[0].size(); j++)
			os << "\t" << tva_[i][j] / (2 * j + 1);
		os << endl;
	}

	os << endl << "Поток энергии в зависимости отположения" << endl << endl;
	for (i = -nr_; i < nr_; i++)
		os << "\t" << (i + 0.5) * dr_;
	os << endl;
	for (i = 0; i < (int)wvxy_.size(); i++)
	{
		os << (i + 0.5 - nr_) * dr_;
		for (j = 0; j < (int)wvxy_[0].size(); j++)
			os << "\t" << wvxy_[i][j];
		os << endl;
	}

	os << endl << "Поток энергии в зависимости от угла" << endl << endl;
	for (i = -nthet_; i < nthet_; i++)
		os << "\t" << (i + 0.5) * dthet_ * 180 / PI;
	os << endl;
	for (i = 0; i < (int)wvvxvy_.size(); i++)
	{
		os << (i + 0.5 - nthet_) * dthet_ * 180 / PI;
		for (j = 0; j < (int)wvvxvy_[0].size(); j++)
			os << "\t" << wvvxvy_[i][j];
		os << endl;
	}
}
