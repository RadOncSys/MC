// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcphysics.h"

class mcRng;

class mcPhysicsCharged : public mcPhysics
{
public:
	mcPhysicsCharged(void);
	virtual ~mcPhysicsCharged(void);

protected:
	// Physics utilities
	static void MoliereScatter(mcRng& rng
		, double scaledLittleB
		, double chi_cc
		, geomVector3D& u
		, double eTotal
		, double betaSquared
		, double effPathLength);

	static double ReducedPathLength(double stepNext, double tscat);
};
