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
	/// –ассто€ние от указанной точки до цилиндра в указанном единичным вектором направлении
	/// до бесконечного цилиндра радиуса R, ось которого совпадает с осью Z.
	/// </summary>
	static double getDistanceToInfiniteCylinderInside(const geomVector3D& p, const geomVector3D& v, double r);
	static double getDistanceToInfiniteCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r);

	/// <summary>
	/// ¬ отличие от предыдущего бесконечного цилиндра этот имеет торцы, 
	/// определ€емые плоскост€ми Z=0 и Z=h
	/// </summary>
	static double getDistanceToCylinderInside(const geomVector3D& p, const geomVector3D& v, double r, double h);
	static double getDistanceToCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r, double h);

	/// <summary>
	/// –ассто€ние от указанной точки до паралелепипеда со сторонами ax и ay и высотой h,
	/// зажатого между плоскост€м Z = 0 и Z = h.
	/// </summary>
	static double getDistanceToPrismInside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h);
	static double getDistanceToPrismOutside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h);

	/// <summary>
	/// ѕересечение с круглым конусом
	/// ‘окус должен быть положительным, т.е. конус должен быть ориентирован так,
	/// что его остра€ часть должна указывать в положительном направлении Z.
	/// </summary>
	static double getDistanceToConeInside(const geomVector3D& p, const geomVector3D& v, double r, double f);
	static double getDistanceToConeOutside(const geomVector3D& p, const geomVector3D& v, double r, double f);

	/// <summary>
	/// ѕересечение с конусообразной боковой поверхностью между z1 b z2 и r1 и r2.
	/// ¬ажно: должно выполн€тьс€ условие z1 < p.z < z2.
	/// </summary>
	static double getDistanceToConeSlabInside(const geomVector3D& p, const geomVector3D& v, double z1, double z2, double r1, double r2);
	static double getDistanceToConeSlabOutside(const geomVector3D& p, const geomVector3D& v, double z1, double z2, double r1, double r2);

	/// <summary>
	/// ѕересечение с ассиметричной бесконечной плоскопараллельной трубой
	/// </summary>
	static double getDistanceToRectanglePipeInside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double y1, double y2);
	static double getDistanceToRectanglePipeOutside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double y1, double y2);

	/// <summary>
	/// ѕересечение с конусом с пр€моугольным симметричным сечением
	/// </summary>
	static double getDistanceToRectangleConeInside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double cosx, double sinx, 
		double y1, double y2, double cosy, double siny, double z);
	static double getDistanceToRectangleConeOutside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double cosx, double sinx,
		double y1, double y2, double cosy, double siny, double z);

	/// <summary>
	/// –асчет длины трека в вокселе.
	/// ¬ аргументах сначала идут координаты начала и конца трека,
	/// затем координаты вершина параллелепипеда.
	/// </summary>
	static double TrLenInVoxel(double x1, double y1, double z1
		, double x2, double y2, double z2
		, double xv1, double yv1, double zv1
		, double xv2, double yv2, double zv2);

	/// <summary>
	/// –ассто€ние от указанной точки до сферы в указанном единичным вектором направлении
	/// до сферы радиуса R, ось которого совпадает с осью Z.
	/// </summary>
	static double getDistanceToSphereInside(const geomVector3D& p, const geomVector3D& v, double r);
	static double getDistanceToSphereOutside(const geomVector3D& p, const geomVector3D& v, double r);

	/// <summary>
	/// –ассто€ние до выпуклого объекта образованном вращением полигона
	/// </summary>
	static double getDistanceToConvexPolygonCircleInside(const geomVector3D& p, const geomVector3D& v, const std::vector<double>& pz, const std::vector<double>& pr);
	static double getDistanceToConvexPolygonCircleOutside(const geomVector3D& p, const geomVector3D& v, const std::vector<double>& pz, const std::vector<double>& pr);
};


