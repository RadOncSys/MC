// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"
#include <vector>

// ����� MLC.
// ��� ������ ��� ���������.
// ������ ������������� ���� �������������� �����
class mcTransportMLC : public mcTransport
{
public:
	mcTransportMLC(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		double focus, double r, double h, double w, double l, 
		double fsx1, double fsx2, double fsy1, double fsy2);
	virtual ~mcTransportMLC(void);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	// ����� ���������� � ����������� �� ����� ���������� � ����������� ������� 
	// ���������� � ������� �����, ������������� � �������� ����������� ����� � ������ X=0
	// ��� ������������ ������ ����� �����.
	double getDistanceToBlockInside(const geomVector3D& p, const geomVector3D& u, double y1, double y2, double xx) const;

	// ���������� ���������� �� MLC ��� ���������� �� ������� ����������
	double getDistanceToMlcOutsideField(const geomVector3D& p, const geomVector3D& u) const;

	// ���������� ���������� �� MLC ��� ���������� ������ ������������� ����
	double getDistanceToMlcInsideField(const geomVector3D& p, const geomVector3D& u) const;

	// ���������� ���������� �� MLC ��� ���������� ������ �������� ������ ����� ������������� ����.
	// ����� �� ���������� �� �������������, ���� ������� �������� � ����� ��������� ����� ���
	// ���������� ������� ����������� (XZ). ����� �� ��������� ������� ������ ���������� ���������� �������.
	double getDistanceToMlcInsideFieldBlock(const geomVector3D& p, const geomVector3D& u, double x1, double x2, double y1, double y2) const;

	// ���������� ���������� �� ������, ������������� ������� ��� ���������� ����������� ������.
	// ���� ����������� � ������� ������������� �� ��������� ���� ���, �� ������������ �������������.
	double getDistanceToMlcRectangleConeInside(const geomVector3D& p, const geomVector3D& u, double x1, double x2, double y1, double y2) const;

	// �������� ����� �� ����� �� ����������� ����������� � ��������� Z = const
	bool isInZPlaneHit(double x, double y) const;

	// ��������������, ��� ������� ��������� ������ ���� �����������.
	// ������������, ��������� �� ��� �� ��������� �������� ���������.
	bool isOutsideMLC(double x, double y) const;

	// ��������������, ��� ����� �� ����� �� ������ XZ � � ���� MLC (�� Z).
	// ����������, ����� �� ����� �� ��������.
	bool isOnOutsideLeaf(const geomVector3D& p) const;

	// ����������, ����� �� ����� �� ��������.
	// ������ ����� ��������, ��� ��� �������������, ��� ��� ��� ���������� � ����� �������.
	bool isOnOnTheLeftLeaf(const geomVector3D& p, double y0, double xx) const;

	void dumpVRMLMLC(ostream& os) const;

	double focus_;	// ���������� �� ��������� �� ������� ���������� �����������
	double r_;		// ������ �������� ����� ���������
	double h_;		// ������ ��������
	double w_;		// ������ ����� ���������
	double l_;		// ����� ���������� ��������
	double fsx1_;	// ������ ����
	double fsx2_;	// (���������� ���������� �������� ����������� ����� ��������,
	double fsy1_;	// �.�. ��� ��������� ��� ����������� ������� ������� ����
	double fsy2_;	// �� ������ ���������)

	// ��������������� ���������� ��� ���������
	double r2_;
	double dx_;		// ���������� ����������� ����� �� �������� ����������� �� X
	double sfy_;	// ���������� ����������� �������� �� Y � �������� �� ������ �����������

	geomVector3D pymin_, pymax_;	// ��������� ������� (0, +/-w/2, 0)
	geomVector3D py1_, py2_;		// ��������� ������� (0, +fsy1/fsy2, 0)
	geomVector3D zxn1_, zxn2_;		// ������� � ������ XZ
};
