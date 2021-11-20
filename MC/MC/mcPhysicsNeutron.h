// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2021] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcPhysicsCharged.h"

class mcPhysicsNeutron : public mcPhysics
{
public:
	mcPhysicsNeutron(void);
	virtual ~mcPhysicsNeutron(void);

	// ������� ���� �� ����������� �������������� (�������� ������� ���������, ���������� ���������).
	// ����������� ������ ������� ���� �� ������.
	double MeanFreePath(double ke, const mcMedium& med, double dens) const override;

	double DoInterruction(mcParticle* p, const mcMedium* med) const override;

	// �������� ������� ������� � ������������ � de/dx � ����������� �������� ������� �� ������� step.
	double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const override;

	bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const override;
};
