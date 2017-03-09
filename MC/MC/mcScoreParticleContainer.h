// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс сбора статистики потока излучения в плоскости.
// Регистрируются координаты точки пересечения плоскости, радиус в плоскости XY и угол по отношению к оси Z.
class mcScoreParticleContainer : public mcScore
{
public:
	mcScoreParticleContainer(const char* module_name, int nThreads);
	virtual ~mcScoreParticleContainer(void);

	void ScoreFluence(const mcParticle& particle) override;

	void setParticleTypeFilter(mc_particle_t t) { ptypeFilter_ = t; }
	double etotal_other() const;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreParticleContainer&);

protected:
	mc_particle_t ptypeFilter_;
	double* etotal_other_;
	struct PlaneParticleRecord
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

	std::vector<std::vector<PlaneParticleRecord>> particles_;
};
