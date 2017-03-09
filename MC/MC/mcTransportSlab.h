// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта в бесконечном по XY слое 
class mcTransportSlab : public mcTransport
{
public:
	mcTransportSlab();
	mcTransportSlab(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double h);
	virtual ~mcTransportSlab(void);

	void setGeometry(double h) { h_ = h; }
	double getHeight() const { return h_; }

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	double h_;
};
