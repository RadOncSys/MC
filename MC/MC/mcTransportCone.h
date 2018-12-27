// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"
#include "../geometry/vec2d.h"

// ����� ���������� � ���������� ��������, 
// ������������ ������������ ��� � ������� ����������� �������
class mcTransportCone :public mcTransport
{
public:
	mcTransportCone();
	mcTransportCone(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r, double h);
	virtual ~mcTransportCone(void);

	void setGeometry(double r, double h);
	double getRadius() const { return r_; }
	double getHeight() const { return h_; }

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDNearInside(const geomVector3D& p) const override;

protected:
	double r_;
	double h_;
	// ����� ����� ��������
	double s_;
	// ������� �������� ���� ������� ������
	double A_;
	// ��������� ������ � ���������� ����� ��������
	geomVector2D g_;
};
