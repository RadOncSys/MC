// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// �������� ���������� � ����������� ������ � ������� �������� ������,
// ���������� ���������� ����� ��������� ������������� ���������� �� �������.
// �������� �������� ������������ ���� ������� ���������� ������ �� ������.
class mcSourceXraySigmaR : public mcSource
{
public:
	mcSourceXraySigmaR(const char* name, int nThreads, mc_particle_t type, double ke, double z, double r, double theta);
	virtual ~mcSourceXraySigmaR(void);

	void sample(mcParticle& p, mcThread* thread) override;
	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceXraySigmaR& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.z_ << endl;
		os << "ROTATION ANGLE = \t" << s.theta_ << endl;
		return os;
	}

protected:
	mc_particle_t type_;
	int q_;			// �����
	double ke_;		// �������
	double z_;		// ��������� ��������� ������� �����
	double sigma_;	// ����� �������� �� ����� ���
	double theta_;	// ���� ������� ���������� �� ������ (� ��������), ������������� ���������

	// ��������������� ����������
	double tr_;		// ���� ������� � ��������
	double uz_;		// �������� ������� �������� �� ��� z
	double sinu_;	// ����� ���� ������� ������� ����������
};
