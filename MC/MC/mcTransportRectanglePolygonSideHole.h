// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2019] Gennady Gorlachev (ggorlachev@roiss.ru) 
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
		double dx, double dy, const std::vector<double>& z, const std::vector<double>& x, const std::vector<double>& y);
	virtual ~mcTransportRectanglePolygonSideHole(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	void SetFieldSize(double x1, double x2, double y1, double y2);

protected:
	void dumpVRMLPolygonSideHole(ostream& os) const;

	// Расчет расстояний в направлении движения, когда частица находится внутри слоя объекта
	double getDistanceOutsideWithinLayer(const geomVector3D& p, const geomVector3D& v) const;

protected:
	// Толщины камней (расстояние от внутренней грани до внешней), коллимирующих поле.
	double dx_;
	double dy_;

	// Координаты полигона грани шторок
	std::vector<double> z_;		// Положения точек полигона по оси Z
	std::vector<double> x_;		// Положения точек полигона шторок X
	std::vector<double> y_;		// Положения точек полигона шторок Y

	// Раскрытие шторок для обеспечения ассиметричного поля.
	// Это физические смешения полигонов внутренней части коллиматора.
	// Положение задается извне в физических см 
	// (сам коллиматор ничего не знает о поле на уровне изоцентра)
	double fsx1_;
	double fsx2_;
	double fsy1_;
	double fsy2_;

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
