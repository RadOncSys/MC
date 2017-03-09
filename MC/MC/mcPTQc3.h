// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"

// ����� ���������� � ������� QC (Standard Imaging), ������������ ��� ������������ ���������� ������.
class mcPTQc3 : public mcTransport
{
public:
	mcPTQc3();
	mcPTQc3(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x,
		int nc, double dr, double ds, double h);
	virtual ~mcPTQc3(void);

	void setGeometry(int nc, double dr, double ds, double h);
	int Nc() const { return nc_; }
	double Dr() const { return dr_; }
	double Ds() const { return ds_; }
	double getHeight() const { return h_; }

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	// ������ ����������� ��� �������������� ������ ���������� �������,
	// ��������� � ���������� �����.

	// ���������� �����
	int nc_;

	// ������� ������
	double dr_;

	// ��� ���������� ����� �� �������
	// ���������� ������� ������� ������ ��������� �� ���������� ds_ �� ������;
	// �������� 2*ds_ � �.�.
	// � ������ ��������� ���������� ������� ������� dr_.
	double ds_;

	// ������ �����
	double h_;
};
