// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2021] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcPhysicsCharged.h"

class mcMediumProton;
class mcRng;

class mcPhysicsProton : public mcPhysicsCharged
{
public:
	mcPhysicsProton(void);
	virtual ~mcPhysicsProton(void);

	// Средний путь до дискретного взаимодействия (рождение второго электрона, тормозного излучения).
	// Непрерывные потери энергии сюда не входят.
	double MeanFreePath(double ke, const mcMedium& med, double dens) const override;

	double DoInterruction(mcParticle* p, const mcMedium* med) const override;

	// Изменяет энергию частицы в соответствии с de/dx и расчитывает передачу энергии на отрезке step.
	double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const override;

	bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const override;

protected:

	static void createnewparticlewithEA(mcRng& rng, mcParticle* p, const mcMediumProton* pmed, int endfID, int pID);
};
