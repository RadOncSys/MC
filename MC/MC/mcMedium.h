// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcElement.h"
#include <string>
#include <vector>
using namespace std;

class mcMedium
{
public:
	mcMedium(void);
	virtual ~mcMedium(void);

	enum STATUS { EMPTY, LOADED, FAILED };

	virtual void read(istream& is) = 0;

protected:
	// Функция извлекает 2 слова из строки.
	// Первое слово сравнивается с именем параметра.
	// При несовпадении устанавливается исключение.
	// По окончании сторка line содержит остаток строки.
	string ParseLine(string& line, const char* param);

public:
	// General:
	STATUS status_;     // Текущее состояние данных
	string name_;       // Name of medium
	double density_;    // Default density of medium
	vector<mcElement> elements_;
};
