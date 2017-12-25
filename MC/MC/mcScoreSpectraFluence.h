// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс скоринга мпектров потока излучения в кольцах.
// Двуменрная матрица по энергии и радиусу накапливает частицы по их количеству (умноженному на вес).
// Параллельно поддерживается матрица того же размера, но накапливающая сразу энергетические спектры
// (чтобы исключить артефакты дискретизации по нергии).
// Параллельно собирается поток энергии в тех же радиальных бинах.
// Учитываются только частицы указанного типа и имеющие кинетическую энергию выше указанного порога.
// При выводе статистики счетчики делятся на площадь кольца в см^2.
class mcScoreSpectraFluence : public mcScore
{
public:
	mcScoreSpectraFluence(const char* module_name, int nThreads, mc_particle_t pt, double ecut, int nr, int ne, double rmax, double emax);
	virtual ~mcScoreSpectraFluence();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreSpectraFluence&);

protected:
	double eFluenceBin(int ir) const;
	double nSpectrumBin(int ir, int ie) const;
	double eSpectrumBin(int ir, int ie) const;

	mc_particle_t pt_;
	int nr_;
	double rstep_;
	double rmax_;
	int ne_;
	double estep_;
	double emax_;
	double ecut_;	// энергия, ниже которой частицы исключаются из статистики
	std::vector<std::vector<double>> energy_fluence_;
	std::vector<std::vector<std::vector<double>>> number_spectra_;
	std::vector<std::vector<std::vector<double>>> energy_spectra_;
};
