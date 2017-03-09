// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс сбора статистики потока фотонного излучения на поверхности сферы.
// Регистрируются координаты точки выхода и угол по отношению к направлению из центра сферы в точку выхода.
class mcScoreSphereFluence : public mcScore
{
public:
	mcScoreSphereFluence(const char* module_name, int nThreads);
	virtual ~mcScoreSphereFluence();

	void ScoreFluence(const mcParticle& particle) override;
	
	double etotal_other() const;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreSphereFluence&);

protected:
	double* etotal_other_;
	struct SFParticleRecord
	{
		float x;
		float y;
		float z;
		float ux;
		float uy;
		float uz;
		float e;
		float r;
		float a;
	};
	
	std::vector<std::vector<SFParticleRecord>> particles_;
};
