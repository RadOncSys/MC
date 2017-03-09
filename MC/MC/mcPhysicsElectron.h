// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcPhysicsCharged.h"

class mcMediumXE;

class mcPhysicsElectron : public mcPhysicsCharged
{
public:
	mcPhysicsElectron(void);
	virtual ~mcPhysicsElectron(void);

	// Средний путь до дискретного взаимодействия (рождение второго электрона, тормозного излучения).
	// Непрерывные потери энергии сюда не входят.
	double MeanFreePath(double ke, const mcMedium& med, double dens) const override;

	double DoInterruction(mcParticle* p, const mcMedium* med) const override;

	// Изменяет энергию частицы в соответствии с de/dx и расчитывает передачу энергии на отрезке step.
	double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const override;

	bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const override;

protected:
	static void MollerScatter(mcParticle* p, const mcMediumXE* med);
	static void DoBremsstrahlung(mcParticle* p, const mcMediumXE* med);

	static void GetMollerEnergies(mcRng& rng, const mcMediumXE* med, double ke_in, double* ke_elec1, double* ke_elec2);
	static double BremsPhotonEnergy(mcRng& rng, const mcMediumXE* med, double ke_in);
	static double BremsPhotonAngle(mcRng& rng, const mcMediumXE* med, double ke_in, double e_phot);
	static void SetUpBremsPhoton(mcParticle* pPhoton, const mcMediumXE* med, double ke_in, double e_phot);
	static double BremsRejectionFunction(double x, double zFactor, double r, double val1, double val2, double val3);
};