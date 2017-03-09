// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"
#include <iostream>
#include "../geometry/vec2d.h"

using namespace std;

//����� ����������� � ������������ �������� ��������
//������� ������� � ���� ������� ���������� �������� (z, r)
class mcTransportPolygon :
	public mcTransport
{
public:
	mcTransportPolygon(void);
	//np - ����� �����(��������), �������� �������
	mcTransportPolygon(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, int np, double h, double r);
	virtual ~mcTransportPolygon(void);

	void setGeometry(int np, double h, double r);

	//virtual void dump(ostream& os) const; 
	//virtual void dumpVRML(ostream& os)const;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	geomVector2D* _ptr;
	int _np;//���������� �����,�������� �������
	double _r;//������ ��������� ��������
	double _h;//������ ��������
};

