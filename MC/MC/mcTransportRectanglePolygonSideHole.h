// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2019] Gennady Gorlachev (ggorlachev@roiss.ru) 
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
		double dx, double dy, const std::vector<double>& z, const std::vector<double>& x, const std::vector<double>& y);
	virtual ~mcTransportRectanglePolygonSideHole(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	void SetFieldSize(double x1, double x2, double y1, double y2);

protected:
	void dumpVRMLPolygonSideHole(ostream& os) const;

	// ������ ���������� � ����������� ��������, ����� ������� ��������� ������ ���� �������
	double getDistanceOutsideWithinLayer(const geomVector3D& p, const geomVector3D& v) const;

protected:
	// ������� ������ (���������� �� ���������� ����� �� �������), ������������� ����.
	double dx_;
	double dy_;

	// ���������� �������� ����� ������
	std::vector<double> z_;		// ��������� ����� �������� �� ��� Z
	std::vector<double> x_;		// ��������� ����� �������� ������ X
	std::vector<double> y_;		// ��������� ����� �������� ������ Y

	// ��������� ������ ��� ����������� �������������� ����.
	// ��� ���������� �������� ��������� ���������� ����� �����������.
	// ��������� �������� ����� � ���������� �� 
	// (��� ���������� ������ �� ����� � ���� �� ������ ���������)
	double fsx1_;
	double fsx2_;
	double fsy1_;
	double fsy2_;

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
