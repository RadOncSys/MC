// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта в шаре с учетом возможности вложения объектов
class mcETransportSphere : public mcTransport
{
public:
	mcETransportSphere();
	mcETransportSphere(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r);

	void setRadius(double r) { r_ = r; }
	double getRadius() const { return r_; }

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

protected:
	double r_;
};
