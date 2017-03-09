// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "../geometry/vec3d.h"

//Purpose:  Базовый класс вспомогательных классов для расчета 
//расстояний до граней, ограничивающих объемы.
class mcGeomSide
{
public:
	virtual ~mcGeomSide() {}

	virtual double getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const = 0;
	virtual double getDNear(const geomVector3D& p) const = 0;
	virtual void dump(ostream& os) const = 0;
};
