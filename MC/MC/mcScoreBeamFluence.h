// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include "mcScoreBeamFluenceSubsource.h"

class mcScoreBeamFluence : public mcScore
{
public:
	mcScoreBeamFluence(const char* module_name, int nThreads, int nr, double rmax);
	virtual ~mcScoreBeamFluence();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreBeamFluence&);

	void addSubsource(mcScoreBeamFluenceSubsource*);

protected:
	int nsub_;
	mcScoreBeamFluenceSubsource* subsources_[10];

	// Нужно для vrml
	int     nr_;
	double  rmax_;
	double  rstep_;
};
