// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

//  ласс транспорта в цилиндрическом кольце, грани которого параллельны
// и внутренн€€ направлена в указанный фокус.
// ¬ собственной системе координат центр последней 
// находитс€ в плоскости с отверстием большего радиуса.
// ‘окус находитс€ на оси Z в положительной стороне.
class mcTransportConicalRing : public mcTransport
{
public:
	mcTransportConicalRing();
	mcTransportConicalRing(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r0, double r1, double h, double f);
	virtual ~mcTransportConicalRing(void);

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
	double f_;  // фокусное рассто€ние внутреннего конуса
	double f_ext_;  // фокусное рассто€ние внешнего конуса конуса

	// —лужебные переменные
	double cosr_;
};
