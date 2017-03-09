// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "../geometry/mtrx3d.h"
#include "mcGeomSide.h"

//Purpose:  Базовый класс геометрических объектов с произвольной формой,
//задаваемой списком ограничивающих его граней.
class mcGeomShape
{
public:
	mcGeomShape();
	virtual ~mcGeomShape();

	virtual double getDistanceInside(const geomVector3D& p, const geomVector3D& v) const;
	virtual double getDistanceOutside(const geomVector3D& p, const geomVector3D& v) const;
	virtual double getDNearInside(const geomVector3D& p) const;

	friend ostream& operator << (ostream& os, const mcGeomShape& g)
	{
		os << "Number of sides:\t" << g.nSides_ << endl;
		for (unsigned int i = 0; i < g.nSides_; i++) g.Sides_[i]->dump(os);
		return os;
	}

protected:
	// Инициализация матрицы преобразования координат.
	// Пользователь может переопределить преобразование координат.
	// Рекомендуемый путь - задание точки начала координат и векторов
	// осей X и Y системы объекта в мировой системе с последующим
	// вызовом данной функции инициализации базового класса.
	virtual void initTransformation() = 0;

	inline double getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const;
	inline double getDNear(const geomVector3D& p) const;

protected:
	// Массив граней, формирующих тело клина
	unsigned int nSides_;
	mcGeomSide** Sides_;
};
