// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// ����� ������������� ����� ������ ��������������.
// � ������ �����������, ��������, ������������� ���������� ������.
// � ����������� ������� ��������� ����� ��������� ��������� � ��������� ������� ���������.
class mcTransportRectangleRing : public mcTransport
{
public:
	mcTransportRectangleRing(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		double x1, double x2, double y1, double y2, double d, double h);
	virtual ~mcTransportRectangleRing(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	bool isXYInside(double x, double y) const;
	bool isXYInsideField(double x, double y) const;
	bool isXYOutsideCollimator(double x, double y) const;
	double distanceToRectangleInside(double x, double y, double X1, double X2, double Y1, double Y2) const;
	double distanceToRectangleOutside(double x, double y, double X1, double X2, double Y1, double Y2) const;

protected:
	double d_;  // ������� ������ (������ ������ �������)
	double h_;  // ������ �������
	double x1_; // ��������� ������ X1
	double x2_; // ��������� ������ X2
	double y1_; // ��������� ������ Y1
	double y2_; // ��������� ������ Y2
};
