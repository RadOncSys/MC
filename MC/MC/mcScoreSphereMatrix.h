// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс измерения дозового распределения в сферическом кольце.
// Изначально предназначен для расчета радиационной защиты головок радиотерапевтических аппаратов.
// Детектор имеет только один слой по радиусу, но по поверхности это двумерная матрица
// (т.е. предусмотрена цилиндрическая несимметричность из-за прямоугольных шторок.
// Деление на ячейки как на глобусе по широтам и меридианам.
// Учитывая сокращение площадей к полюсам (ось Z) разбиение на широты нелинейное 
// подобно переменному радиусу в радиальных кольцах на плоскости.
class mcScoreSphereMatrix : public mcScore
{
public:
	mcScoreSphereMatrix(const char* module_name, int nThreads, int np, int nm, double rmin, double rmax);
	virtual ~mcScoreSphereMatrix();

	void ScoreFluence(const mcParticle& particle) override;
	void ScorePoint(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0) override;

	void ScoreLine(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0
		, const geomVector3D& p1) override;

	virtual void scoreEnergyInVoxel(int iThread, int ip, int im, double edep);

	void CE2D() override;

	double	Dose(int ip, int im) const;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreSphereMatrix&);

protected:
	// Расстояние из точки p0 в направлении v до конуса, соответствующего широте angle
	double distanceToParallel(const geomVector3D& p0, const geomVector3D& v, int ip, bool isPOutside);

	// Расстояние из точки p0 в направлении v до плоскости, соответствующей долготе angle
	double distanceToMeridian(const geomVector3D& p0, const geomVector3D& v, int im, bool isMPositive);

protected:
	int np_;				// количество параллелей
	int nm_;				// количество меридианов
	double rmin_, rmax_;	// радиусы сфер, между которыми слой детектора
	double* MAll_;
	double** M_;

	double r2min_, r2max_;
	double pstep_;
	double mstep_;
};
