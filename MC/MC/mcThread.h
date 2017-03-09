// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcRng.h"

class mcParticle;

class mcThread					//what 4 this class
{
public:
	mcThread();
	~mcThread(void);

	int id() const { return threadIdx_; }
	void setId(int id);

	mcRng& rng() { return rng_; }

	mcParticle* NextParticle();
	void RemoveParticle();
	mcParticle* DuplicateParticle();

	mcParticle** CurrentParticle() { return &pCurParticle_; }
	mcParticle* ParticleStack() { return particleStack_; }

protected:
	int threadIdx_;
	mcRng rng_;
	mcParticle* pCurParticle_;
	mcParticle* particleStack_;
};
