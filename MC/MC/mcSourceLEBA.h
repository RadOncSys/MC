// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// ���� ����� �������������� �������
enum spectrum_distr_t { SPECTRUM_GAUSS = 0, SPECTRUM_TRIANGLE, SPECTRUM_PRISM };

// ���� ���������� ��������
enum profile_distr_t { PROFILE_PRISM = 0, PROFILE_GAUSS, PROFILE_EXPONENT };

// �������� ���������� �� ���������� �� ������� ��������� ������ �������.
// �� ���������� ���������������� �������������� ������������� �� �������
// � ������� �������������� �������� � ����� �������������� ��� ���������� �� ���� �������������.
// � ������ ��������� �������� ������ ���� �����������.
// ������� ������� ������ ������������ �� �� ����� ��������.
// LEBA - Linear Electron Beam Accelerator.
class mcSourceLEBA : public mcSource
{
public:
	// name - ������� ������;
	// sptype - ��� ����� ������� (0 - �������, 1 - �����������, 2 - ������������);
	// kemean - ������� ������� � ��� (�������� �������������� ������� ��� ��� ������ ����);
	// kewidth -������ ��������������� ������� � ��� (������ ��������������, ������ ��������� ������������, ��� 2 ����� ���������);
	// z - ��������� �������� � ��;
	// rsigma - ���������� ���������� ������������� ������������� �� �������.
	mcSourceLEBA(const char* name, int nThreads, spectrum_distr_t sptype, profile_distr_t prf_type,
		double kemean, double kewidth, double z, double rsigma);

	void sample(mcParticle& p, mcThread* thread) override;
	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceLEBA& s)
	{
		os << (const mcSource&)s;
		os << "POSITION = \t" << s.z_ << endl;
		os << "TYPE = \t" << s.sptype_ << endl;
		os << "KEMEAN = \t" << s.kemean_ << endl;
		os << "KEWIDTH = \t" << s.kewidth_ << endl;
		os << "RSIGMA = \t" << s.rsigma_ << endl;
		return os;
	}

protected:
	spectrum_distr_t sptype_;
	profile_distr_t prf_type_;
	int q_;				// �����
	double kemean_;		// ������� �������
	double kewidth_;	// ������ ��������������� �������
	double z_;			// ��������� ��������� ������� �����
	double rsigma_;		// ���������� ���������� ������������� �� �������
};
