// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

class mcPhaseSpaceIO;

// Класс регистрации частиц в файле фазового пространства
class mcScorePHSP : public mcScore
{
public:
	mcScorePHSP(const char* module_name, const char* fname);
	virtual ~mcScorePHSP(void);

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

protected:
	mcPhaseSpaceIO* phsp_;
};
