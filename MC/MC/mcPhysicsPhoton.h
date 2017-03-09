// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcphysics.h"

class mcMediumXE;
class mcRng;

class mcPhysicsPhoton : public mcPhysics
{
public:
	mcPhysicsPhoton(void);
	virtual ~mcPhysicsPhoton(void);

	double MeanFreePath(double ke, const mcMedium& med, double dens) const override;
	double DoInterruction(mcParticle* p, const mcMedium* med) const override;
	double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const override;
	bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const override;

protected:
	static void ProducePair(mcParticle* p, const mcMediumXE* pMedium);
	static void GetPairEnergies(mcRng& rng, const mcMediumXE* pMedium, double e_in, double* ke_elec1, double* ke_elec2);
	static void GetPairAngle(mcRng& rng, const mcMediumXE* pMedium, double eScaled_in, double ke_elec, double* cosTheta, double* sinTheta);
	static double PairRejectionFunction(double x, double zFactor, double r, double oneOverR, double val1);

	static void ComptonScatter(mcParticle* p, const mcMediumXE* pMedium);
	static void GetComptonSplit(mcRng& rng, double e_in,
		double* e_phot, double* cosTheta_phot, double* sinTheta_phot,
		double* ke_elec, double* cosTheta_elec, double* sinTheta_elec);

	static double DoPhotoelectric(mcParticle* p, const mcMediumXE* pMedium);
	static void GetPhotoelectricEnergies(double rnum, const mcMediumXE* pMedium, double e_in, double* e_dep, double* ke_elec, double* e_fluor);
	static void ChangePhotoelectronDirection(mcRng& rng, double ke_elec, geomVector3D& u);

	static void RayleighScatter(mcParticle* p, const mcMediumXE* pMedium);
	static void GetRayleighAngle(mcRng& rng, const mcMediumXE* pMedium, double e_phot, double* cosTheta, double* sinTheta);
};
