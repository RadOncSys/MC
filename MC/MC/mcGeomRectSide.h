// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcgeomside.h"

//Purpose:  ��������������� ����� ��� ������� ���������� ��.
//��������������, ����������� �������������� � ������������ � 
//��������������� ����� ����������� �������.
class mcGeomRectSide : public mcGeomSide
{
public:
	mcGeomRectSide(const geomVector3D& p
		, const geomVector3D& Vx
		, const geomVector3D& Vy
		, double ax, double ay);

	double getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const override;
	double getDNear(const geomVector3D& p) const override;
	void dump(ostream& os) const override {}

protected:
	// ��������� ���������� �����
	geomVector3D P_;      // ���������� �������� �������
	geomVector3D Vx_;     // ������� ������ �� ������� (� ������� �������)
	geomVector3D Vy_;
	double ax_;         // ������� ������ ��������������
	double ay_;

	// ����������� ���������
	geomVector3D N_;      // ������� � ��������� ��������������
};
