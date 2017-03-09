// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта в цилиндрическом кольце
class mcTransportRing : public mcTransport
{
public:
	mcTransportRing();
	mcTransportRing(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r0, double r1, double h);
	virtual ~mcTransportRing(void);

	void setGeometry(double r0, double r1, double h);
	double R0() const { return r0_; }
	double R1() const { return r1_; }
	double getHeight() const { return h_; }

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
};
