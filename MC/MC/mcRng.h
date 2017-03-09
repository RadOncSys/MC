// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcRng
{
public:
	mcRng(void);
	mcRng(int ijSeed, int klSeed);
	~mcRng(void);

	void init(int ijSeed, int klSeed);

	double rnd();

protected:
	double  genArray_[97];
	double* p96Gen_;   // Pointer to last element of generator array
	double* pIGen_;
	double* pJGen_;

	double decrement_;
};
