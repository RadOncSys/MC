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

	double kEmax(void)const { return (double)dedx1_proto.size(); }
	virtual void read(istream& is);
	void createDB();
	const double AtomicWeight() const;	// ������� ��� �����, �/����

private:
	//--------------------------------
	// ��������� ������ (������!)
	//--------------------------------
	double gdEdxStragglingGaussVarianceConstPart();	// ���������� � ���������� ���������� (�� ������� � ����) ����� �������� ��������� ����������� �������� dE/dx. 
	double gRadiationLength();	// ���������� � ���������� �������� �������� ������������ ����� �������� �������� rpp-2006-book.pdf 27.4.1 p.263 (eq.27.23) ��� ������� ������������ ����� ���������� �������� �������� InverseRadiationLength  (�������� ����������� (� ������ 2007 ���� ��� ����������� Dahl'�))
	void gSigmaInelastic(int Ap = 1, int Zp = 1);	// ���������� �������� ������� ���������� ��������������

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

	std::shared_ptr<std::vector<std::shared_ptr<mcEndfP>>> ENDFdata;
	double microsigmaforelement(int A, int Z, double kE) const;
};
