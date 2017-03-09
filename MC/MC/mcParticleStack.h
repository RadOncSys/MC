// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcRng;
class mcParticle;

class mcParticleStack
{
public:
	mcParticleStack(int nThreads = 1);
	~mcParticleStack(void);

protected:
	int nThreads_;
	mcRng* pRng_;
	mcParticle** pCurParticle_;
	mcParticle** particleStack_;
};
