// Radiation Oncology Monte Carlo open source project
//
// Author: [2020] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"

// Класс транспорта в горизонтально движущейся пары шторок коллиматора 
// с гранями, постоянно смотрящими в фокус пучка.
// Ориентация граней расчитывается автоматически в конструкторе 
// исходя из размера поля и фокусного расстояния.
// Изначально геометрия создана для расчета системы коллимации с фокусированными парами шторок
// и опробована на описании экранов этажерки коллиматора Terabalt.
class mcTransportJawPairFocused : public mcTransport
{
public:
	mcTransportJawPairFocused();
	// sad - source isocenter distance,
	// scd_ - расстояние отисточника доо основания шторки, 
	// h - толщина шторки, 
	// dx - ширина прямоугольной части (на уровне основания), 
	// dy - ширина в перпендикулярном направлении
	mcTransportJawPairFocused(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, 
		double sad, double scd, double h, double dx, double dy);
	virtual ~mcTransportJawPairFocused(void);

	void setFS(double x1, double x2);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	// Геометрия аппарата
	double sad_;
	double scd_;

	// Толщина шторки
	double h_;

	// Размер шторки на уровне основания
	double dx_;
	double dy_;

	// Размер ассиметричного поля в физических единицах
	double fsx1_, fsx2_;
	double fscenter_;

	// Служебные параметры геометрии для ускорения расчетов
	double x1l_;
	double x2l_;
	double x1r_;
	double x2r_;
	geomVector3D nl_;
	geomVector3D nr_;
	geomVector3D pl_;
	geomVector3D pr_;
};
