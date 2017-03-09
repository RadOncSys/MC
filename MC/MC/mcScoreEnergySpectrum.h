// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс скоринга энергетическог спектра частиц определенного типа.
// Изначально создан скорее для тестовых задач, чем для реальных расчетов.
class mcScoreEnergySpectrum : public mcScore
{
public:
	mcScoreEnergySpectrum(const char* module_name, int nThreads, mc_particle_t pt, int ne, double emax, double rmax);
	virtual ~mcScoreEnergySpectrum();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreEnergySpectrum&);

protected:
	mc_particle_t pt_;
	int ne_;
	double emax_;
	double estep_;
	double rmax_;
	std::vector<std::vector<double>> espec_;
};
