#include "mcScoreBeamFluenceXY.h"
#include "mcParticle.h"
#include "mcThread.h"
#include "mcTransport.h"

mcScoreBeamFluenceXY::mcScoreBeamFluenceXY(const char* module_name, int nThreads,
	int nx, int ny, double psx, double psy, int nebins, double emax)
	:mcScore(module_name, nThreads),
	nx_(nx), ny_(ny), psx_(psx), psy_(psy), nebins_(nebins), emax_(emax)
{
	minx_ = -0.5 * nx_ * psx_;
	miny_ = -0.5 * ny_ * psy_;
	estep_ = emax_ / nebins_;

	int len = nThreads * nx_ * ny_;
	intencity_all_ = new double[len];
	memset(intencity_all_, 0, len * sizeof(double));
	intencity_ = new double*[nThreads];
	intencity_[0] = intencity_all_;
	for (int i = 1; i < nThreads_; i++)
		intencity_[i] = intencity_[i - 1] + nx_ * ny_;

	len = nThreads * nebins_;
	spectrum_all_ = new double[len];
	memset(spectrum_all_, 0, len * sizeof(double));
	spectrum_ = new double*[nThreads];
	spectrum_[0] = spectrum_all_;
	for (int i = 1; i < nThreads_; i++)
		spectrum_[i] = spectrum_[i - 1] + nebins_;
}

mcScoreBeamFluenceXY::~mcScoreBeamFluenceXY()
{
	if (intencity_all_) delete[] intencity_all_;
	if (intencity_) delete[] intencity_;
	if (spectrum_all_) delete[] spectrum_all_;
	if (spectrum_) delete[] spectrum_;
}

void mcScoreBeamFluenceXY::ScoreFluence(const mcParticle& particle)
{
	if (particle.t != MCP_PHOTON)
		return;

	int iThread = particle.thread_->id();
	double edep = particle.ke * particle.weight;
	etotal_[iThread] += edep;

	// !!! Scoring вызывается до перемещения частицы на поверхность объекта к которому привязана частица.
	// Нужно ее перенести здесь на поверхность.
	geomVector3D p = particle.p + (particle.u * (-particle.p.z() / particle.u.z()));

	double dx = p.x() - minx_;
	double dy = p.y() - miny_;
	if (dx < 0 || dy < 0) return;
	int ix = int(dx / psx_);
	int iy = int(dy / psy_);
	if (ix >= nx_ || iy >= ny_) return;
	intencity_[iThread][iy*nx_ + ix] += edep;

	// Спектр
	int idx = (int)(particle.ke / estep_);
	if (idx < nebins_)
		spectrum_[iThread][idx] += edep;
}

double mcScoreBeamFluenceXY::Intencity(int iThread, int ix, int iy) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluenceXY::Intencity: thread exceed container size");
	if (ix >= 0 && ix < nx_ && iy >= 0 && iy < ny_)
		return intencity_[iThread][iy*nx_ + ix];
	else
		return 0;
}

double mcScoreBeamFluenceXY::Intencity(int ix, int iy) const
{
	double f = 0;
	if (ix >= 0 && ix < nx_ && iy >= 0 && iy < ny_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += intencity_[i][iy*nx_ + ix];
	}
	return f;
}

double mcScoreBeamFluenceXY::Spectrum(int iThread, int idx) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluenceXY::Intencity: thread exceed container size");
	if (idx >= 0 && idx < nebins_)
		return spectrum_[iThread][idx];
	else
		return 0;
}

double mcScoreBeamFluenceXY::Spectrum(int idx) const
{
	double f = 0;
	if (idx >= 0 && idx < nebins_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += spectrum_[i][idx];
	}
	return f;
}

void mcScoreBeamFluenceXY::dumpVRML(ostream& os) const
{
	os << "# mcScoreBeamFluenceXY Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	//const geomMatrix3D& mttow = transport_->MT2W();

	//os << "Shape {" << endl;
	//os << "  appearance Appearance {" << endl;
	//os << "    material Material {" << endl;
	//os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	//os << "    }" << endl;
	//os << "  }" << endl;
	//os << "  geometry IndexedLineSet {" << endl;

	//os << "    coord Coordinate {" << endl;
	//os << "      point [" << endl;

	//// Концентрические круги
	//for(ir=1; ir<=m_nr; ir++) {
	//	double r = m_rstep * ir;
	//	for(it=0; it<360; it+=da) {
	//		geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), 0) * mttow;
	//		os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	//		p = geomVector3D(r*sin(mPi*(it+da)), r*cos(mPi*(it+da)), 0) * mttow;
	//		os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	//		count++;
	//	}
	//}

	//os << "      ]" << endl;
	//os << "    }" << endl;

	//os << "    coordIndex [" << endl;
	//for(it=0; it<count; it++)
	//	os << "      " << 2*it << ' ' << 2*it+1 << " -1" << endl;
	//os << "    ]" << endl;
	//os << "  }" << endl;
	//os << "}" << endl;
}

void mcScoreBeamFluenceXY::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	os << endl;
	os << "Fluence intencity Map" << endl;
	os << endl;

	int i, j;
	for (i = 0; i < nx_; i++)
		os << "\t" << minx_ + (i + 0.5) * psx_;
	os << endl;

	for (j = 0; j < ny_; j++)
	{
		os << miny_ + (j + 0.5) * psy_;
		for (i = 0; i < nx_; i++)
			os << "\t" << Intencity(i, j);
		os << endl;
	}
	os << endl;

	os << "Mean spectrum" << endl;
	for (i = 0; i < nebins_; i++)
		os << "\t" << EnergyBin(i);
	os << endl;
	for (i = 0; i < nebins_; i++)
		os << "\t" << Spectrum(i);
	os << endl;
	os << endl;
}
