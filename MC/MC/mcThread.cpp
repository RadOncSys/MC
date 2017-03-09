#include "mcThread.h"
#include "mcParticle.h"

#define MC_STACK_SIZE   15

mcThread::mcThread() : threadIdx_(0)
{
	particleStack_ = new mcParticle[MC_STACK_SIZE];
	pCurParticle_ = particleStack_ - 1;
}

mcThread::~mcThread(void)
{
	delete[] particleStack_;
}

void mcThread::setId(int id)
{
	threadIdx_ = id;
	rng_.init(33, id + 37);
}

mcParticle* mcThread::NextParticle()
{
	pCurParticle_++;
	if (pCurParticle_ >= particleStack_ + MC_STACK_SIZE)
		throw std::exception("mcThread::NextParticle: particle stack overflow");
	return pCurParticle_;
}

void mcThread::RemoveParticle()
{
	if (pCurParticle_ < particleStack_)
		throw std::exception("mcThread::RemoveParticle: remove not existing particle");
	pCurParticle_--;
}

mcParticle* mcThread::DuplicateParticle()
{
	if (pCurParticle_ >= particleStack_ + MC_STACK_SIZE - 1)
		throw std::exception("mcThread::DuplicateParticle: can not duplicate particle, stack will oferflow");
	if (pCurParticle_ < particleStack_)
		throw std::exception("mcThread::DuplicateParticle: can not duplicate particle, nothing to duplicate");
	mcParticle* p = pCurParticle_;
	pCurParticle_++;

	// Просто удобное место для сброса MFPs.
	// Это нужно сделать, так как размножение частиц всешда взаимодействие, после которого все начинается сначала
	p->mfps = 0;
	*pCurParticle_ = *p;

	return pCurParticle_;
}
