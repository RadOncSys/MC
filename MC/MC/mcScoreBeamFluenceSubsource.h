// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcParticle.h"

// Служебный класс для реализации в одном score возможность отслеживания подисточников в одной симуляции.
// В этом классе содержится информация для фильтрации частиц по типу и положению последнего взаимодействия.
// Здесь же накапливается вся информация о потоке от подисточника.
class mcScoreBeamFluenceSubsource
{
public:
	mcScoreBeamFluenceSubsource(const char* name, int nThreads,
		enum mc_particle_t, double zmin, double zmax, double focus,
		int nr, double rmax,
		int ne, double emax,
		int na, double amax);
	~mcScoreBeamFluenceSubsource(void);

	void ScoreFluence(const mcParticle& particle, const geomVector3D& p);
	void dumpStatistic(ostream&) const;
	friend ostream& operator << (ostream&, const mcScoreBeamFluenceSubsource&);

public:
	int    Nr() const { return nr_; }
	double MaxR() const { return rmax_; }
	double StepR() const { return rstep_; }

	int    Ne() const { return ne_; }
	double MaxE() const { return emax_; }
	double StepE() const { return estep_; }

	int    Na() const { return na_; }
	double MaxA() const { return amax_; }
	double StepA() const { return astep_; }

	double Intencity(int iThread, int ir) const;
	double Intencity(int ir) const;

	double Spectrum(int iThread, int ir, int ie) const;
	double Spectrum(int ir, int ie) const;
	double SpectrumSigma(int ir, int ie) const;

	double Angle(int iThread, int ir, int ia) const;
	double Angle(int ir, int ia) const;
	double AngleSigma(int ir, int ia) const;

protected:
	char name_[32];
	enum mc_particle_t ptype_;
	double zmin_;
	double zmax_;
	int nThreads_;
	double focus_;

	int     nr_;
	double  rmax_;
	double  rstep_;
	double** intencity_;
	double* intencity_all_;

	int     ne_;
	double  emax_;
	double  estep_;
	double** spectra_;
	double* spectra_all_;

	int     na_;
	double  amax_;
	double  astep_;
	double** angle_;
	double* angle_all_;
};
