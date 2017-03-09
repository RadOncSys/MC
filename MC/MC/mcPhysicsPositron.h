// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcPhysicsElectron.h"

class mcPhysicsPositron : public mcPhysicsElectron
{
public:
	mcPhysicsPositron(void);
	virtual ~mcPhysicsPositron(void);

	double MeanFreePath(double ke, const mcMedium& med, double dens) const override;
	double DoInterruction(mcParticle* p, const mcMedium* med) const override;

	// Изменяет энергию частицы в соответствии с de/dx и расчитывает передачу энергии на отрезке step.
	double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const override;

	bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const override;

protected:
	static void BhabhaScatter(mcParticle* p, const mcMediumXE* med);
	static void AnnihilateInFlight(mcParticle* p);
	static void AnnihilateAtRest(mcParticle* p);
	static void GetBhabhaEnergies(mcRng& rng, const mcMediumXE* pMedium, double ke_in, double* ke_posi, double* ke_nega);
	static void GetAnnihilationEnergies(mcRng& rng, double ke_in, double p_in, double* e_phot1, double* e_phot2);

};
