// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

// Класс сбора статистики потока излучения в координатах XY.
class mcScoreBeamFluenceXY : public mcScore
{
public:
	mcScoreBeamFluenceXY(const char* module_name, int nThreads,
		int nx, int ny, double psx, double psy,
		int nebins, double emax);
	virtual ~mcScoreBeamFluenceXY();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreBeamFluenceXY&);

public:
	double Intencity(int iThread, int ix, int iy) const;
	double Intencity(int ix, int iy) const;
	double Spectrum(int iThread, int idx) const;
	double Spectrum(int idx) const;
	double EnergyBin(int idx) const { return emax_ * (idx + 0.5) / nebins_; }
	int NEnergyBins() const { return nebins_; }

protected:
	int    nx_;
	int    ny_;
	double psx_;
	double psy_;
	double* intencity_all_;
	double** intencity_;

	// Spectrum
	int nebins_;
	double emax_;
	double estep_;
	double* spectrum_all_;
	double** spectrum_;

	// Служебные переменные
	double minx_, miny_;
};
