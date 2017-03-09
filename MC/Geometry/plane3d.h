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

	// Положение секущей плоскости вдоль оси перепендикулярной плоскости.
	// Изначально предназначено для обслуживания генератора 3D матриц
	// из стэка плоских изображений.
	double getPlanePosition() const;

	const geomMatrix3D& getRtoPMatrix() const { return mRtoP_; }
	const geomMatrix3D& getPtoRMatrix() const { return mPtoR_; }

	// Utilities

	// Возвращает точку пересечения линии, проходящей через указанные точки.
	// Результат в мироаой системе координат.
	// Если линия параллельеа плоскрсти, то вывешивается Exception.
	geomVector3D crossByLine(const geomVector3D& p0, const geomVector3D& p1)const;

	// Определяет точку пересечения отрезка с плоскостью.
	// Возвращает false, если пересечения нет, или оно не между точками.
	// Точка пересечения задается в координатах плоскости.
	bool crossByEdge(const geomVector3D& p0, const geomVector3D& p1, double& x, double& y)const;

	// Ближайшее расстояние от точки до плоскости
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
	geomMatrix3D mRtoP_;	// Преобразование координат из плоскости в мировую систему
	geomMatrix3D mPtoR_;	// Преобразование координат из мировую системs в плоскость
};
