// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcParticle;
class mcMedium;
class geomVector3D;

// ������� ����� ������������� ���������� ���������
class mcPhysics
{
public:
	mcPhysics(void);
	virtual ~mcPhysics(void);

	// ������ ������� ����� ������� � �����.
	virtual double MeanFreePath(double ke, const mcMedium& med, double dens) const = 0;

	// ������������ ���� ��������������.
	// ���������������, ��� ��������� �� ������ �������� ���������� �� ������� ������� (�����).
	// ���� � ���������� �������������� ������ ����� ������� ���������� ���������, ������
	// ���������� � ���� �� ��������� ���������� ������ �������� �������.
	// ������������ ������������ � ����� � ���������� �������������� �������.
	virtual double DoInterruction(mcParticle* p, const mcMedium* med) const = 0;

	// ����������� ������� �� �������� ����������. ���������� ���������� �� ���� �������.
	// ��� ������� ��� ������ ��������� ���������� �������.
	// � ���������� step ������������ �������� ���. ��� ���������� ������ �� ����� ���������� �� ���������,
	// ���� ��� ������������� ����������� ������ �������� ������.
	// � ��������� ������ ��������� ����������, ��� ���������� ������� ����������� �� �������.
	virtual double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const = 0;

	// �������� ������� �� ������� ������� ���� ����������� � ��� ������������� �� �����������.
	// ������������ ������� ���������� � edep.
	// ����������� ������ �����, ��� ��� �������� �������� ���� �����������, ��� �������� ������ ������.
	virtual bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const = 0;

	// Physics utilities
	static void GetRandomPhi(double rnum, double* cosPhi, double* sinPhi);
	static void ChangeDirection(double cosTheta, double sinTheta, double cosPhi, double sinPhi, geomVector3D& u);
	static void GoInRandomDirection(double rnum1, double rnum2, geomVector3D& u);

	static mcParticle* DuplicateParticle(mcParticle* p);
	static void DiscardParticle(mcParticle* p);
};
