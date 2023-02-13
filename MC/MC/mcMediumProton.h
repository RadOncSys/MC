// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2021] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcMedium.h"
#include "mcDefs.h"

// ����� �������� ���������� ���������� ����� ��� ���������� ��������.
class mcMediumProton : public mcMedium
{
public:
  mcMediumProton(void);
  virtual ~mcMediumProton(void);

  const double kEmax(void)const{return (double)dedx1_proto.size();}
  virtual void read(istream& is);
	const double AtomicWeight() const;	// ������� ��� �����, �/����

private:
//--------------------------------
// ��������� ������ (������!)
//--------------------------------
	const double	gdEdxStragglingGaussVarianceConstPart(); // ���������� � ���������� ���������� (�� ������� � ����) ����� �������� ��������� ����������� �������� dE/dx. 
	const double	gRadiationLength();	// ���������� � ���������� �������� �������� ������������ ����� �������� �������� rpp-2006-book.pdf 27.4.1 p.263 (eq.27.23) ��� ������� ������������ ����� ���������� �������� �������� InverseRadiationLength  (�������� ����������� (� ������ 2007 ���� ��� ����������� Dahl'�))
	const void		gSigmaInelastic(int Ap=1, int Zp=1);	// ���������� �������� ������� ���������� ��������������

public:
  // TODO !!!
  // ���� ����������� � ���������� ������ � ��������� PEGS
  // �� ������ ��������� ������������ � ������ ��������.
  // ��� ������� ��������� ������ � �������� �������� ������
  // ������ ����� mcPhysicsProton.

	//vector<double> dEdx_; // dE/dx - ����������� �� �����
	//vector<double> sigma_in_;	// ������� ��������� ������� �������������� ��� ������� 1,2...,
	
	// �� ��������� �� ������� � ���� ����� ��������� �������� (sigma^2) dE/dx
	double dEdxStragglingGaussVarianceConstPart_;
	
	
  //// �������� ��������, ��� � EGS � �������� ������� ������ � ��������������� ��������.
  //// ������ � ����� �� �� �� ��� ��������, ���� ������� ��������� dE/dX � �������� ��������?

  //// Log energy index:
  //double iLogKE0_proto;
  //double iLogKE1_proto;

  //// Transport:
  double radLength;          // Radiation length, [cm](!)
  //double pStep;

  //// Moliere multiple scattering:
  //double scaledLittleB;
  //double chi_cc;

  //// Interactions:
  //// � ������ ���������� ������������ Branching ����������
  //// ������������� ����������� ���� ������������ ���������� �������.
  //vector<double> br10_proto;         // Branching ratio 1
  //vector<double> br11_proto;

  //// Transport:
  // ������������ ������������ S(E)=S0[Ei]+S1[Ei]*E, Ei=int(E)
  // � ���������� ������ ������� ������������ ��������
  // �� ���� ������� �������� ������������
  vector<double> sigma0_proto;       // Fictitious cross section, 1/cm
  vector<double> sigma1_proto;
  vector<double> dedx0_proto;        // Linear energy loss rate (stopping power)(Mev/cm)
  vector<double> dedx1_proto;
  //vector<double> stepSize0_proto;    // Step size
  //vector<double> stepSize1_proto;

  double transCutoff_proto; // Energy cutoff for proton transport
};
