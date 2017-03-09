// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcgeomside.h"
#include "../geometry/vec2d.h"

//Purpose:  ¬спомогательный класс дл€ расчета рассто€ни€ до.
//треугольника, произвольно расположенного в пространстве и 
//представл€ющего грань трехмерного объекта.
class mcGeomTriangleSide : public mcGeomSide
{
public:
	mcGeomTriangleSide(const geomVector3D& p
		, const geomVector3D& Vx
		, const geomVector3D& Vy
		, double ax, double ay);

	double getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const override;
	double getDNear(const geomVector3D& p) const override;

	void dump(ostream& os) const override {}

protected:
	// ѕараметры сообщаемые извне
	geomVector3D P_;        // координаты основной вершины
	geomVector3D Vx_;       // векторы сторон из вершины (в мировой системе)
	geomVector3D Vy_;
	double ax_;             // размеры сторон треугольника
	double ay_;

	// ѕроизводные параметры
	geomVector3D N_;        // нормаль к плоскости треугольника
	geomVector3D VRy_;      // ось Y системы треугольника в мировой системе

	geomVector2D p_[3];     // вершин треугольника
	geomVector2D v_[3];     // векторы ребер
	geomVector2D n_[3];     // левые нормали ребер

	double a_[3];           // длины ребер
};
