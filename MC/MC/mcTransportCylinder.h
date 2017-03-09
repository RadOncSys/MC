// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта в цилиндре
class mcTransportCylinder : public mcTransport
{
public:
	mcTransportCylinder();
	mcTransportCylinder(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r, double h);
	virtual ~mcTransportCylinder(void);

	void setGeometry(double r, double h);
	double getRadius() const { return r_; }
	double getHeight() const { return h_; }

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

	//protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	double r_;
	double h_;
};
