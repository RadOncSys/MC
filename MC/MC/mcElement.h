// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcElement
{
public:
	mcElement(void) : atomicNumber(0), atomicMass(0), partsByNumber(0)
	{
		atomicSymbol[0] = static_cast<char>(0);
	}
	~mcElement(void) {}

	char atomicSymbol[3];
	int atomicNumber;
	double atomicMass;
	double partsByNumber;
};
