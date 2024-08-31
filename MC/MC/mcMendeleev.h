// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <vector>

// Класс меток о том для каких элементов таблицы менделеева нужна данные и для каких загружены
class mcMendeleev
{
public:
	mcMendeleev();
	mcMendeleev(const mcMendeleev& t);

	std::vector<bool> IsNecessary;
	std::vector<bool> IsLoad;
};
