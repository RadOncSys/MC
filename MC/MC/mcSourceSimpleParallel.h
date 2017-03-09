// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// ������������ (z=const) ������������ (v) �������� ������ (type), 
// � �������� ������������� �� ������ (��� sigmaE=0 - ����������������)
// ������������� ������� ��������������� � 0.
// ��� �������� ��������� ������������� rx=ry
class mcSourceSimpleParallel : public mcSource
{
public:
	mcSourceSimpleParallel(void);
	mcSourceSimpleParallel(const char* name, int nThreads,
		mc_particle_t type, double ke, const geomVector3D& p, const geomVector3D& v,
		double rx, double ry, double sigmaKE = 0);
	virtual ~mcSourceSimpleParallel(void);

	void init(mc_particle_t type        // ��� ������
		, double ke                // ������������ �������, ��� (�������-!)
		, const geomVector3D& p    // ����� �������� ������
		, const geomVector3D& v	// ����������� �������� (��������� ������)
		, double rx				// x-������ ������ ���������
		, double ry				// y-������ ������ ���������
		, double sigmaKE			// ����� ��������� �������
	);
	void sample(mcParticle& p, mcThread* thread) override;

	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceSimpleParallel& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.p_ << endl;
		os << "DIRECTION = \t" << s.v_ << endl;
		os
			<< "RX = \t" << s.rx_ << endl
			<< "RY = \t" << s.ry_ << endl
			<< "SIGMA_KE = \t" << s.sigmaKE_ << endl;
		return os;
	}

protected:
	mc_particle_t type_;
	double ke_;
	geomVector3D p_;
	geomVector3D v_;
	int q_; // �����
	double rx_;				// x-������ ������ ���������
	double ry_;				// y-������ ������ ���������
	double sigmaKE_;			// ����� ��������� �������
};
