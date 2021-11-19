// Radiation Oncology Monte Carlo open source project
//
// Author: [2020] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"

// ����� ���������� � ������������� ���������� ���� ������ ����������� 
// � �������, ��������� ���������� � ����� �����.
// ���������� ������ ������������� ������������� � ������������ 
// ������ �� ������� ���� � ��������� ����������.
// ���������� ��������� ������� ��� ������� ������� ���������� � ��������������� ������ ������
// � ���������� �� �������� ������� �������� ����������� Terabalt.
class mcTransportJawPairFocused : public mcTransport
{
public:
	mcTransportJawPairFocused();
	// sad - source isocenter distance,
	// scd_ - ���������� ����������� ��� ��������� ������, 
	// h - ������� ������, 
	// dx - ������ ������������� ����� (�� ������ ���������), 
	// dy - ������ � ���������������� �����������
	mcTransportJawPairFocused(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, 
		double sad, double scd, double h, double dx, double dy);
	virtual ~mcTransportJawPairFocused(void);

	void setFS(double x1, double x2);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	// ��������� ��������
	double sad_;
	double scd_;

	// ������� ������
	double h_;

	// ������ ������ �� ������ ���������
	double dx_;
	double dy_;

	// ������ �������������� ���� � ���������� ��������
	double fsx1_, fsx2_;
	double fscenter_;

	// ��������� ��������� ��������� ��� ��������� ��������
	double x1l_;
	double x2l_;
	double x1r_;
	double x2r_;
	geomVector3D nl_;
	geomVector3D nr_;
	geomVector3D pl_;
	geomVector3D pr_;
};
