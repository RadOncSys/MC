// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"

// ����� ���������� � ������������� ���������� ���� ������ ����������� � �������� �������� ������.
// ������ ����� ������������ �����.
// ��� X ������ � ���������� �� ��������.
// ���������� ��������� ������� ��� ������� ������� ���������� ��������� ������ ����������.
class mcTransportJawPairRounded : public mcTransport
{
public:
	mcTransportJawPairRounded();
	// r - ������ ����������� �����, h - ������� ������, dx - ������ ������������� �����
	mcTransportJawPairRounded(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r, double h, double dx);
	virtual ~mcTransportJawPairRounded(void);

	void setFS(double x1, double x2) { fsx1_ = x1; fsx2_ = x2; fscenter_ = 0.5 * (x1 + x2); }

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;
	
	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceToLeftOutside(const geomVector3D& p, const geomVector3D& u) const;

	void dumpVRMLJaw(ostream& os, bool isLeft) const;

protected:
	// ������ ����������� �����
	double r_, r2_;

	// ������� ������
	double h_;

	// ������ ������������� �����
	double dx_;

	// ������ �������������� ���� � ���������� ��������
	double fsx1_, fsx2_;
	double fscenter_;
};
