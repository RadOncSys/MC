// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "vec3d.h"

class geomPlane3D
{
public:
	geomPlane3D();
	geomPlane3D(const geomPlane3D&);
	geomPlane3D(const geomVector3D& p, const geomVector3D& n, const geomVector3D& xv) { set(p, n, xv); }
	void set(const geomVector3D& p, const geomVector3D& n, const geomVector3D& xv);

	const geomVector3D&	getPoint() const { return p_; }
	const geomVector3D&	getNormal() const { return n_; }
	const geomVector3D&	getXAxis() const { return xv_; }
	bool isAxial() const { return n_.x() == 0 && n_.y() == 0; }

	void setPoint(const geomVector3D& p) { p_ = p; }
	void setNormal(const geomVector3D& n) { n_ = n; }
	void setXAxis(const geomVector3D& xv) { xv_ = xv; }
	double posZ() const { return p_.z(); }

	// ��������� ������� ��������� ����� ��� ����������������� ���������.
	// ���������� ������������� ��� ������������ ���������� 3D ������
	// �� ����� ������� �����������.
	double getPlanePosition() const;

	const geomMatrix3D& getRtoPMatrix() const { return mRtoP_; }
	const geomMatrix3D& getPtoRMatrix() const { return mPtoR_; }

	// Utilities

	// ���������� ����� ����������� �����, ���������� ����� ��������� �����.
	// ��������� � ������� ������� ���������.
	// ���� ����� ����������� ���������, �� ������������ Exception.
	geomVector3D crossByLine(const geomVector3D& p0, const geomVector3D& p1)const;

	// ���������� ����� ����������� ������� � ����������.
	// ���������� false, ���� ����������� ���, ��� ��� �� ����� �������.
	// ����� ����������� �������� � ����������� ���������.
	bool crossByEdge(const geomVector3D& p0, const geomVector3D& p1, double& x, double& y)const;

	// ��������� ���������� �� ����� �� ���������
	double nearestDistance(const geomVector3D& p)const;

	friend istream& operator >> (istream&, geomPlane3D&);
	friend ostream& operator << (ostream&, const geomPlane3D&);

protected:
	void initMRtoP();

protected:
	geomVector3D p_;		// point, belonging the plane
							// usualy the origion of the plane coordinates
	geomVector3D n_;		// normal vector to the plane
	geomVector3D xv_;		// Vector, that specifies x-axis of the plane
	geomMatrix3D mRtoP_;	// �������������� ��������� �� ��������� � ������� �������
	geomMatrix3D mPtoR_;	// �������������� ��������� �� ������� ������s � ���������
};
