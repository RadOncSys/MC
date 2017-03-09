// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"
#include <vector>

// ����� �������������� ��������� � ���������������.
// ����� ������ ��� �������� ��������� ����������� ������-�.
// � ����������� ������� ��������� ����� ��������� ��������� � ��������� �������� ���������.
// ��� Z ���������� �� ��������� � ������� �������� (��� ������� ���������).
// ��������� ������ ����������� ��� ������������.
// ������������� ������ ����� �������������� ��������� ������ ������� ��������� �������.
class mcTransportRectanglePolygonSideHole : public mcTransport
{
public:
	mcTransportRectanglePolygonSideHole(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		double dx, double dy, std::vector<double>& z, std::vector<double>& x, std::vector<double>& y);
	virtual ~mcTransportRectanglePolygonSideHole(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	void dumpVRMLPolygonSideHole(ostream& os) const;

protected:
	double dx_; // ��������� ������� ������ ������ X ������������ ������ ��������� �������
	double dy_; // ��������� ������� ������ ������ Y ������������ ������ ��������� �������

	// ���������� �������� ����� ������
	std::vector<double> z_;		// ��������� ����� �������� �� ��� Z
	std::vector<double> x_;		// ��������� ����� �������� ������ X
	std::vector<double> y_;		// ��������� ����� �������� ������ Y

	// ��������� ���������
	int nlayers_;				// ���������� ����� �������
	double xmax_;
	double ymax_;
	double zmax_;				// ������� ����������� ����� ��� Z
	std::vector<double> cosx_;	// �������� ���������� ������, ����������� ��� Y ������������ ��� Z
	std::vector<double> sinx_;	// ������ (������������� ��� ������������� �����)
	std::vector<double> cosy_;	// �������� ���������� ������, ����������� ��� X ������������ ��� Z
	std::vector<double> siny_;	// ������ (������������� ��� ������������� �����)
};
