// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс скоринга энергетическог спектра частиц определенного типа.
// Изначально создан скорее для тестовых задач, чем для реальных расчетов.
class mcScoreEnergyFluence : public mcScore
{
public:
	mcScoreEnergyFluence(const char* module_name, int nThreads, mc_particle_t pt, int nr, double rmax);
	virtual ~mcScoreEnergyFluence();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreEnergyFluence&);

protected:
	mc_particle_t pt_;
	int nr_;
	double rstep_;
	double rmax_;
	std::vector<std::vector<double>> fluence_;
};
