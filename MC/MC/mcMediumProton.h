// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcEndfP.h"
#include "mcMedium.h"
#include "mcDefs.h"

// ����� �������� ���������� ���������� ����� ��� ���������� ��������.
class mcMediumProton : public mcMedium
{
public:
	mcMediumProton(void);
	virtual ~mcMediumProton(void);

	const double kEmax(void)const { return (double)dedx1_proto.size(); }
	virtual void read(istream& is);
	void createDB();
	const double AtomicWeight() const;	// ������� ��� �����, �/����

private:
	//--------------------------------
	// ��������� ������ (������!)
	//--------------------------------
	const double gdEdxStragglingGaussVarianceConstPart();	// ���������� � ���������� ���������� (�� ������� � ����) ����� �������� ��������� ����������� �������� dE/dx. 
	const void	 gSigmaInelastic(int Ap = 1, int Zp = 1);	// ���������� �������� ������� ���������� ��������������
	const double gRadiationLength();

public:
	// �� ��������� �� ������� � ���� ����� ��������� �������� (sigma^2) dE/dx
	double dEdxStragglingGaussVarianceConstPart_;
	
	//// Transport:
	double radLength;          // Radiation length, [cm](!)

	// Transport:
	// ������������ ������������ S(E)=S0[Ei]+S1[Ei]*E, Ei=int(E)
	// � ���������� ������ ������� ������������ ��������
	// �� ���� ������� �������� ������������
	vector<double> sigma0_proto;	// Fictitious cross section, 1/cm (���� � ������� ���������������)
	vector<double> sigma1_proto;
	vector<double> dedx0_proto;		// Linear energy loss rate (stopping power)(Mev/cm)
	vector<double> dedx1_proto;

	double transCutoff_proto;		// Energy cutoff for proton transport

	vector<mcEndfP> ENDFdata;
	double microsigmaforelement(int A, int Z, double kE) const;
};
