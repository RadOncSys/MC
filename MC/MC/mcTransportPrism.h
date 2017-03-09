// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"

// Класс транспорта в параллелепипеде.
// Центр системы координат объекта находится
// в середине грани основания (Z=const)
class mcTransportPrism : public mcTransport
{
public:
	mcTransportPrism(void);
	mcTransportPrism(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		double ax, double ay, double az);
	virtual ~mcTransportPrism(void);

	void setGeometry(double ax, double ay, double az);
	double ax() const { return ax_; }
	double ay() const { return ay_; }
	double az() const { return az_; }

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	double ax_;
	double ay_;
	double az_;
};
