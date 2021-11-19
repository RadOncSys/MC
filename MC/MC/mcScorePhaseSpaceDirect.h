// Radiation Oncology Monte Carlo open source project
//
// Author: [2020] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

/// <summary>
/// Класс сбора потока частиц в радиально симметричной системе.
/// Обслуживает радиально симметричную модель источника, состоящую из набора частиц, 
/// являющегося продуктом прямой симуляции части радиотерапевтического аппарата.
/// </summary>
class mcScorePhaseSpaceDirect : public mcScore
{
public:
	/// <summary>
	/// Конструктор класса сбора статистикти по распределению частиц по энергии положению и углам
	/// в радиальносимметричной системе.
	/// </summary>
	/// <param name="module_name">имя модуля, с которым ассоциирован данный scoring</param>
	/// <param name="nThreads">количество потоков</param>
	/// <param name="ptype">тип частиц, для которых собирается статистика</param>
	/// <param name="emax">максимальная энергия в МэВ</param>
	/// <param name="rmax">максимальный радиус в см</param>
	/// <param name="stat_file">имя файла вывода результатов</param>
	mcScorePhaseSpaceDirect(const char* module_name, int nThreads, enum mc_particle_t ptype, 
		double emax, double rmax, const char* stat_file);
	virtual ~mcScorePhaseSpaceDirect();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScorePhaseSpaceDirect&);

protected:
	enum mc_particle_t ptype_;
	double emax_;
	double rmax_;
	std::string model_file_;

	double energyScale_;
	double rScale_;				// наиболее неопределенная привязка
	double radialAngleScale_;	// компоненты Ux векторов скорости ([-1,1] -> [0,2])
	double azimutAngleScale_;	// компоненты Uy векторов скорости ([-1,1] -> [0,2])

	// Массивы накопления частиц (три измерения - threads, particles, parameters)
	unsigned short* data_;
	std::vector<unsigned> threadIndexes_;

	// Debug
	unsigned nr_;
	double dr_;
	std::vector<double> weight_;
	std::vector<double> fluence_;
	std::vector<double> rangle_;
	std::vector<double> aangle_;
	std::vector<double> rangle2_;
	std::vector<double> aangle2_;
};
