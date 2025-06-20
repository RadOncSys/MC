// Radiation Oncology Monte Carlo open source project
//
// Author: [2025] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

// Класс скоринга нейтральных частиц для поддержки парадигмы расчета распределений потоков
// Матрица 4D. Четвертая  координата - логарифм энергии.
class mcScoreFluence3D : public mcScore
{
public:
	/// <summary>
	/// Сетка 3D матрицы задается как симметричная по всем осям.
	/// Координата z0 как обычно указывает привязку центра матрицы к объекту транспорта,
	/// например глубина залегания центра матрицы в водном фантоме.
	/// Следует помнить о размерах, так как легко вывалиться за пределы оперативной памяти.
	/// Класс предназначен для расчета потоков нейтронов и фотонов.
	/// Разбиение по энергиям внутри на данном этапе определно статически.
	/// </summary>
	mcScoreFluence3D(const char* module_name, int nThreads, 
		int nx, int ny, int nz, double psx, double psy, double psz, double z0);
	virtual ~mcScoreFluence3D(void);

	// Скоринг потока ведется в единицах длины треков в ячейки для каждого энергетического интервала.
	// Подобно дозам флюенс расчитывается в постобработке.
	void ScoreLine(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0
		, const geomVector3D& p1) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

protected:
	int    nx_, ny_, nz_;
	double z0_;
	double psx_;
	double psy_;
	double psz_;
	double** M_;
};
