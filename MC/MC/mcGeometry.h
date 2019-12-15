// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <vector>

class geomVector3D;

class mcGeometry
{
public:
	/// <summary>
	/// ���������� �� ��������� ����� �� �������� � ��������� ��������� �������� �����������
	/// �� ������������ �������� ������� R, ��� �������� ��������� � ���� Z.
	/// </summary>
	static double getDistanceToInfiniteCylinderInside(const geomVector3D& p, const geomVector3D& v, double r);
	static double getDistanceToInfiniteCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r);

	/// <summary>
	/// � ������� �� ����������� ������������ �������� ���� ����� �����, 
	/// ������������ ����������� Z=0 � Z=h
	/// </summary>
	static double getDistanceToCylinderInside(const geomVector3D& p, const geomVector3D& v, double r, double h);
	static double getDistanceToCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r, double h);

	/// <summary>
	/// ���������� �� ��������� ����� �� �������������� �� ��������� ax � ay � ������� h,
	/// �������� ����� ���������� Z = 0 � Z = h.
	/// </summary>
	static double getDistanceToPrismInside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h);
	static double getDistanceToPrismOutside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h);

	/// <summary>
	/// ����������� � ������� �������
	/// ����� ������ ���� �������������, �.�. ����� ������ ���� ������������ ���,
	/// ��� ��� ������ ����� ������ ��������� � ������������� ����������� Z.
	/// </summary>
	static double getDistanceToConeInside(const geomVector3D& p, const geomVector3D& v, double r, double f);
	static double getDistanceToConeOutside(const geomVector3D& p, const geomVector3D& v, double r, double f);

	/// <summary>
	/// ����������� � �������������� ������� ������������ ����� z1 b z2 � r1 � r2.
	/// �����: ������ ����������� ������� z1 < p.z < z2.
	/// </summary>
	static double getDistanceToConeSlabInside(const geomVector3D& p, const geomVector3D& v, double z1, double z2, double r1, double r2);
	static double getDistanceToConeSlabOutside(const geomVector3D& p, const geomVector3D& v, double z1, double z2, double r1, double r2);

	/// <summary>
	/// ����������� � ������������� ����������� ������������������ ������
	/// </summary>
	static double getDistanceToRectanglePipeInside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double y1, double y2);
	static double getDistanceToRectanglePipeOutside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double y1, double y2);

	/// <summary>
	/// ����������� � ������� � ������������� ������������ ��������
	/// </summary>
	static double getDistanceToRectangleConeInside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double cosx, double sinx, 
		double y1, double y2, double cosy, double siny, double z);
	static double getDistanceToRectangleConeOutside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double cosx, double sinx,
		double y1, double y2, double cosy, double siny, double z);

	/// <summary>
	/// ������ ����� ����� � �������.
	/// � ���������� ������� ���� ���������� ������ � ����� �����,
	/// ����� ���������� ������� ���������������.
	/// </summary>
	static double TrLenInVoxel(double x1, double y1, double z1
		, double x2, double y2, double z2
		, double xv1, double yv1, double zv1
		, double xv2, double yv2, double zv2);

	/// <summary>
	/// ���������� �� ��������� ����� �� ����� � ��������� ��������� �������� �����������
	/// �� ����� ������� R, ��� �������� ��������� � ���� Z.
	/// </summary>
	static double getDistanceToSphereInside(const geomVector3D& p, const geomVector3D& v, double r);
	static double getDistanceToSphereOutside(const geomVector3D& p, const geomVector3D& v, double r);

	/// <summary>
	/// ���������� �� ��������� ������� ������������ ��������� ��������
	/// </summary>
	static double getDistanceToConvexPolygonCircleInside(const geomVector3D& p, const geomVector3D& v, const std::vector<double>& pz, const std::vector<double>& pr);
	static double getDistanceToConvexPolygonCircleOutside(const geomVector3D& p, const geomVector3D& v, const std::vector<double>& pz, const std::vector<double>& pr);
};


