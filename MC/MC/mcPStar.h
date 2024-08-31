// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
// Classes to manage cross sections load for proton nuclear interractions
// from ENDF format
//---------------------------------------------------------------------------
#pragma once
#include <string>
#include <vector>
#include <memory>
#include "mcMendeleev.h"

// Строка параметров для одной энергии файла PSTAR
struct mcPStarRecord
{
	// T    STOP(e)    STOP(n)    STOP(t)
	double D[4] = { 0, 0, 0, 0 };
};

// Полная таблица тормозных способностей для одного элемента
class mcPStarTable
{
public:
	std::string ElementName;
	std::string InternaName;
	int Z;
	int A;
	std::vector<mcPStarRecord> Data;

	// Извлечение данных для одной энергии с помощью look-up и 
	// линейной интерполяции энергий в логанифмическом масштабе
	mcPStarRecord GetDataForEnergy(double e);

	void LoadFromStream(std::istream& is);
};

// Class that keeps PSTAR elements stopping power database for protons
class mcPStar
{
public:
	mcPStar();
	~mcPStar();

	std::vector<std::shared_ptr<mcPStarTable>> Tables;
	
	std::shared_ptr<mcPStarTable> GetTableForElement(int Z) const;

	void LoadFromPath(const std::string& path, mcMendeleev& table);
};
