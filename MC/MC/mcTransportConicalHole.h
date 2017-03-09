// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс конического отверстия в цилиндре.
// Класс создан после и очень похож на mcTransportConicalRing.
// В собственной системе координат центр последней находится в плоскости с отверстием большего радиуса.
// Фокус находится на оси Z в положительной стороне.
// Поэтому, в типичной конфигурации расчетов ось Z объекта направлена вдоль -Z.
class mcTransportConicalHole : public mcTransport
{
public:
	mcTransportConicalHole();
	mcTransportConicalHole(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r0, double r1, double h, double f);
	virtual ~mcTransportConicalHole(void);

	void setGeometry(double r0, double r1, double h, double f);
	double R0() const { return r0_; }
	double R1() const { return r1_; }
	double getHeight() const { return h_; }
	double F() const { return f_; }

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	double r0_; // внутренний радиус
	double r1_; // внешний радиус
	double h_;  // высота кольца
	double f_;  // фокусное расстояние

	// Служебные переменные
	double cosr0_;
};
