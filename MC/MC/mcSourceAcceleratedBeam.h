// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"
#include <vector>

// �������� � ���� ������������ ����� �� ������ ����������.
// ���� ������ ������ ������������� �������� ���������� �������:
// X Y Z Vx Vy Vz E(MeV) Weight
class mcSourceAcceleratedBeam : public mcSource
{
public:
	mcSourceAcceleratedBeam(const char* name, int nThreads, double z);

	void sample(mcParticle& p, mcThread* thread) override;
	void loadData(istream& is);
	void dumpVRML(ostream& os) const override;

protected:
	// �������� ��������� (����������� ������� � ����������� ������� ���������) �� ��� Z
	double z_;
	int nparticles_;
	std::vector<int> idx_;
	std::vector<mcParticle> particles_;
};
