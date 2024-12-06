// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"


// TODO: Класс стартовал с существующего транспорта mcETransportConvexPolygonCircle.
// Отличие только в том, что ось типа радиально симметричной фигуры изгибается.
// Далее нужно заменить функционал mcETransportConvexPolygonCircle 
// на правильный для mcTransportTube3D




// Класс транспорта в изгибающейся полнотелой трубке типа эндостата брахитерапии.
// При симуляции реального эндостато нужно создавать два объекта 
// этого типа один (воздух) вложенный в другой (материал эндостата).
class mcTransportTube3D : public mcTransport
{
public:
	mcTransportTube3D(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, std::vector<double>& pz, std::vector<double>& pr);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	std::vector<double> pz_;
	std::vector<double> pr_;
};
