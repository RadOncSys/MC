// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcscore.h"

/// <summary>
/// Класс сбора статистики потока частиц в виде карт распределений в радиально симметричной системе.
/// </summary>
class mcScorePhaseSpaceConcentrator : public mcScore
{
public:
	/// <summary>
	/// Конструктор класса сбора статистикти по распределению частиц по энергии положению и углам
	/// в радиальносимметричной системе.
	/// </summary>
	/// <param name="module_name">имя модуля, с которым ассоциирован данный scoring</param>
	/// <param name="nThreads">количество потоков</param>
	/// <param name="ptype">тип частиц, для которых собирается статистика</param>
	/// <param name="focus">расстояние до эффективного фокуса, из которого предположительно влетают частицы (например, центра источника)</param>
	/// <param name="isXray">чисто тормоной спектр, или гамма-источник с фиксированными энергиями нерассеянных фотонов</param>
	/// <param name="ne">количество бинов по энергии</param>
	/// <param name="nr">количество бинов по радиусу</param>
	/// <param name="nthet">количество бинов по углу</param>
	/// <param name="naxial">количество бинов по азимуту</param>
	/// <param name="emax">максимальная энергия в МэВ</param>
	/// <param name="rmax">максимальный радиус в см</param>
	/// <param name="thetmax">максимальный угол в градусах</param>
	/// <param name="stat_file">имя файла вывода результатов</param>
	mcScorePhaseSpaceConcentrator(
		const char* module_name, int nThreads, 
		enum mc_particle_t ptype, bool isXray, double focus,
		int ne, int nr, int nthet, int naxial, 
		double emax, double rmax, double thetmax, 
		const char* stat_file);
	virtual ~mcScorePhaseSpaceConcentrator();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScorePhaseSpaceConcentrator&);

protected:
	// Вычисление индекса 4D массива статистики, куда должна быть зарегистрирована частица.
	int particleIdx(const mcParticle& particle) const;

	// Вывод статистики по двум распределениям sude by side
	// (предназначенный в первую очередь для сравнения симулированных, и сгенерированных моделью)
	void dumpDistributions(ostream& os, double* dsim, double* dmodel) const;

protected:
	enum mc_particle_t ptype_;
	int isxray_;
	int ne_;
	int nr_;
	int nthet_;
	int naxial_;

	double focus_;
	double emax_;
	double rmax_;
	double thetmax_;		// (1 - cos(thetmax_abs_ * Pi / 180))
	double thetmax_abs_;	// максимальный угол в градусах

	std::string model_file_;

	// Служебные переменные
	int nne_;
	int nnr_;
	int nnt_;
	double de_;
	double dr_;
	double dthet_;
	double daxial_;
	double* data_;	// 4-х мерный массив данных дифференциальной гистограммы
};
