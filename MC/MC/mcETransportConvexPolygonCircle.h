// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта в полнотелом выпуклом цилиндрически симметричном объекте,
// образованном вращением полигона
class mcETransportConvexPolygonCircle : public mcTransport
{
public:
	mcETransportConvexPolygonCircle(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, std::vector<double>& pz, std::vector<double>& pr);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	std::vector<double> pz_;
	std::vector<double> pr_;
};
