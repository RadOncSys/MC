#include "mcPhysicsNeutron.h"
#include "mcMediumNeutron.h"
#include "mcPhysicsCommon.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcDefs.h"
#include <float.h>

mcPhysicsNeutron::mcPhysicsNeutron(void)
{
}

mcPhysicsNeutron::~mcPhysicsNeutron(void)
{
}

bool mcPhysicsNeutron::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	//if (p->ke <= ((const mcMediumNeutron.transCutoff_proto)
	//{
	//	edep = p->ke;
		DiscardParticle(p);
		return true;
	//}
	//else
	//	return false;
}

double mcPhysicsNeutron::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	throw std::exception("Not implemented");

	return 0;
}

double mcPhysicsNeutron::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	throw std::exception("Not implemented");

	const mcMediumNeutron& m = (const mcMediumNeutron&)med;
	double e_dep = 0;


	return e_dep;
}

double mcPhysicsNeutron::DoInterruction(mcParticle* p, const mcMedium* med) const
{
	throw std::exception("Not implemented");

	// Возвращаем энергию, выделившуюся в точке.
	return 0;
}
