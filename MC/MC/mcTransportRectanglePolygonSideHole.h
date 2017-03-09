// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"
#include <vector>

// Класс прямоугольного отверстия в параллелепипеде.
// Класс создан для описания основного коллиматора Рокуса-Р.
// В собственной системе координат центр последней находится в плоскости входного отверстия.
// Ось Z направлена от источника в сторону изоцента (для лучшего понимания).
// Положение шторок описывается как симметричное.
// Ассиметричные шторки могут симулироваться смещением центра системы координат объекта.
class mcTransportRectanglePolygonSideHole : public mcTransport
{
public:
	mcTransportRectanglePolygonSideHole(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		double dx, double dy, std::vector<double>& z, std::vector<double>& x, std::vector<double>& y);
	virtual ~mcTransportRectanglePolygonSideHole(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	void dumpVRMLPolygonSideHole(ostream& os) const;

protected:
	double dx_; // положение внешних граней шторок X относительно центра координат объекта
	double dy_; // положение внешних граней шторок Y относительно центра координат объекта

	// Координаты полигона грани шторок
	std::vector<double> z_;		// Положения точек полигона по оси Z
	std::vector<double> x_;		// Положения точек полигона шторок X
	std::vector<double> y_;		// Положения точек полигона шторок Y

	// Служебные переменые
	int nlayers_;				// количество слоев объекта
	double xmax_;
	double ymax_;
	double zmax_;				// толщина коллиматора вдоль оси Z
	std::vector<double> cosx_;	// Косинусы плоскостей граней, паралельных оси Y относительно оси Z
	std::vector<double> sinx_;	// Синусы (положительные для расширяющейся части)
	std::vector<double> cosy_;	// Косинусы плоскостей граней, паралельных оси X относительно оси Z
	std::vector<double> siny_;	// Синусы (положительные для расширяющейся части)
};
