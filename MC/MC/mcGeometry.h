// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class geomVector3D;

class mcGeometry
{
public:
	/// <summary>
	/// Расстояние от указанной точки до цилиндра в указанном единичным вектором направлении
	/// до бесконечного цилиндра радиуса R, ось которого совпадает с осью Z.
	/// </summary>
	static double getDistanceToInfiniteCylinderInside(const geomVector3D& p, const geomVector3D& v, double r);
	static double getDistanceToInfiniteCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r);

	/// <summary>
	/// В отличие от предыдущего бесконечного цилиндра этот имеет торцы, 
	/// определяемые плоскостями Z=0 и Z=h
	/// </summary>
	static double getDistanceToCylinderInside(const geomVector3D& p, const geomVector3D& v, double r, double h);
	static double getDistanceToCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r, double h);

	/// <summary>
	/// Расстояние от указанной точки до паралелепипеда со сторонами ax и ay и высотой h,
	/// зажатого между плоскостям Z = 0 и Z = h.
	/// </summary>
	static double getDistanceToPrismInside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h);
	static double getDistanceToPrismOutside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h);

	/// <summary>
	/// Пересечение с круглым конусом
	/// Фокус должен быть положительным, т.е. конус должен быть ориентирован так,
	/// что его острая часть должна указывать в положительном направлении Z.
	/// </summary>
	static double getDistanceToConeInside(const geomVector3D& p, const geomVector3D& v, double r, double f);
	static double getDistanceToConeOutside(const geomVector3D& p, const geomVector3D& v, double r, double f);

	/// <summary>
	/// Пересечение с ассиметричной бесконечной плоскопараллельной трубой
	/// </summary>
	static double getDistanceToRectanglePipeInside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double y1, double y2);
	static double getDistanceToRectanglePipeOutside(const geomVector3D& p, const geomVector3D& v,
		double x1, double x2, double y1, double y2);

	/// <summary>
	/// Пересечение с конусом с прямоугольным симметричным сечением
	/// </summary>
	static double getDistanceToRectangleConeInside(const geomVector3D& p, const geomVector3D& v,
		double x, double cosx, double sinx, double y, double cosy, double siny, double z);
	static double getDistanceToRectangleConeOutside(const geomVector3D& p, const geomVector3D& v,
		double x, double cosx, double sinx, double y, double cosy, double siny, double z);

	/// <summary>
	/// Расчет длины трека в вокселе.
	/// В аргументах сначала идут координаты начала и конца трека,
	/// затем координаты вершина параллелепипеда.
	/// </summary>
	static double TrLenInVoxel(double x1, double y1, double z1
		, double x2, double y2, double z2
		, double xv1, double yv1, double zv1
		, double xv2, double yv2, double zv2);

	/// <summary>
	/// Расстояние от указанной точки до сферы в указанном единичным вектором направлении
	/// до сферы радиуса R, ось которого совпадает с осью Z.
	/// </summary>
	static double getDistanceToSphereInside(const geomVector3D& p, const geomVector3D& v, double r);
	static double getDistanceToSphereOutside(const geomVector3D& p, const geomVector3D& v, double r);
};


