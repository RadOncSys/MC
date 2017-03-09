// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcgeomshape.h"

//Purpose:  ����� �������� ��������� �������� ������������ �����.
//������� ������� ����� - ������������� �����������.
class mcGeomWedgeTriangle : public mcGeomShape
{
public:
	mcGeomWedgeTriangle(double x0, double y0, double x1, double xw);

	friend ostream& operator << (ostream& os, const mcGeomWedgeTriangle& g) {
		os << *(const mcGeomShape*)&g;
		os << "Wede specific:" << endl;
		os << "Polygon (x0,y0,x1):\t" << g.x0_ << '\t' << g.y0_ << '\t' << g.x1_ << endl;
		os << "Wedge length (xw):\t" << g.xw_ << endl;
		return os;
	}

protected:
	double x0_, y0_;    // ���������� ������� �������� �����
	double x1_;         // ��������� ������� �����
	double xw_;         // ����� ����� (� �����������, ���������������� ����������� �����)
};
