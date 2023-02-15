// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
// Classes to manage cross sections load for proton nuclear interractions.
// Current version supports ICRU-63 data.
//---------------------------------------------------------------------------
#pragma once
#include <string>
#include <vector>
#include <memory>

// Ветор спектров родившихся частиц при отдельно взятом угле
class mcCSNuclearForAngleSpectrum
{
public:
	// Угол в градусах, для которого задан данный спектр
	double Angle;

	// Верхние границы диапазона энергий каждого бина спектра
	std::vector<double> SourceBinUps;

	// Интенсивности бинов исходого спектра
	std::vector<double> SourceSpectrum;

	void Clear();
};

// Crossections per istope per incident particle energy
class mcCSNuclearParticleEnergy
{
public:
	// Кинетическая энергия падающего протона
	double Energy;

	// Суммарное сечение ядерных взаимодействий 
	// данной энергии для родительского элемнта (mb)
	double TotalCrossSection;
	double ProtonCrossSection;
	double NeutronCrossSection;

	// Массив спектров рождающихся протонов
	std::vector<mcCSNuclearForAngleSpectrum> ProtonAngles;

	// Массив спектров рождающихся нейтронов
	std::vector<mcCSNuclearForAngleSpectrum> NeutronAngles;

	// Other particles
	// В данной версии не симулируются в предположении локального выделения энергии.
	// Т.е., остаток вероятности суммарного сечения приводит к локальноу 
	// убиванию падающего протон и выделению его кинетической энергии в точке.

	void Clear();
};

// Class that keeps coss sections for one atomic isotope
class mcCSNuclear
{
public:
	// Назавание изотопа, включающее атомное имя и атомный вес.
	// Используется как уникальный идентификатор.
	std::string ElementName;

	// Массивы таблиц сечений отдельно выбранных энергий
	std::vector<mcCSNuclearParticleEnergy> Energies;

	void Clear();
};
