#include "mcParticleStack.h"
#include "mcParticle.h"
#include "mcRng.h"

#define MC_STACK_SIZE   15

mcParticleStack::mcParticleStack(int nThreads) :nThreads_(nThreads)
{
	pRng_ = new mcRng[nThreads];
	pCurParticle_ = new mcParticle*[nThreads];
	particleStack_ = new mcParticle*[nThreads];
	for (int i = 0; i < nThreads; i++)
	{
		pRng_[i].init(33, i + 1);
		particleStack_[i] = new mcParticle[MC_STACK_SIZE];
		pCurParticle_[i] = particleStack_[i];
	}
}

mcParticleStack::~mcParticleStack(void)
{
	delete[] pRng_;
	for (int i = 0; i < nThreads_; i++)
		delete[] particleStack_[i];
	delete[] pCurParticle_;
	delete[] particleStack_;
}
