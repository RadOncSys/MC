// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2021] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcPhysicsCharged.h"
#include <vector>
using namespace std;

class mcMediumProton;
class mcRng;

class mcPhysicsProton : public mcPhysicsCharged
{
public:
	mcPhysicsProton(void);
	virtual ~mcPhysicsProton(void);

	// ������� ���� �� ����������� �������������� (�������� ������� ���������, ���������� ���������).
	// ����������� ������ ������� ���� �� ������.
	double MeanFreePath(double ke, const mcMedium& med, double dens) const override;

	double DoInterruction(mcParticle* p, const mcMedium* med) const override;

	// �������� ������� ������� � ������������ � de/dx � ����������� �������� ������� �� ������� step.
	double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const override;

	bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const override;

protected:

	static void createnewparticleswithEA(mcRng& rng, mcParticle* p, const mcMediumProton* pmed, int endfID, vector<int>* quantity);
	static void getKallbachMannAngle(mcRng& rng, int endfID, mcParticle* p, const mcMediumProton* pmed, int keIN, int eoutID);
};
