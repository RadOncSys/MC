// Radiation Oncology Monte Carlo open source project
//
// Author: [2020] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcETransportTrap.h"

// Transport class, which moves particle to the sphere surface
// and calls fluence scoring if it moves out of the center.
// All particles complete life inside this object.
class mcTransportSphereTrap : public mcETransportTrap
{
public:
	mcTransportSphereTrap(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r);
	~mcTransportSphereTrap(void);

	// Начало транспорта переписано, так как не нужно запускать симуляцию частиц.
	// Их просто нужно отфильтровать.
	void beginTransport(mcParticle& p) override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double r_;	// sphere radius
};
