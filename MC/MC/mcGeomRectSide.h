// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcgeomside.h"

//Purpose:  Вспомогательный класс для расчета расстояния до.
//прямоугольника, произвольно расположенного в пространстве и 
//представляющего грань трехмерного объекта.
class mcGeomRectSide : public mcGeomSide
{
public:
	mcGeomRectSide();
	mcGeomRectSide(const geomVector3D& p,
		const geomVector3D& Vx, const geomVector3D& Vy,
		double ax, double ay);

	void SetGeometry(const geomVector3D& p,
		const geomVector3D& Vx, const geomVector3D& Vy,
		double ax, double ay);

	double getDistance(const geomVector3D& p, const geomVector3D& v) const override;
	double getDNear(const geomVector3D& p) const override;
	void dump(ostream& os) const override {}

protected:
	// Полуразмеры сторон прямоугольника
	double ax_;         
	double ay_;

	// Преобразование координат из внешней системы во внутреннюю
	geomMatrix3D m_; 
};
