// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Ловушка частиц.
// Преобразует координаты в собственную систему и не сдвигая 
// частицу вызывает ассоциированный скоринг.
// Затем частица уничтожается.
class mcETransportTrap : public mcTransport
{
public:
	mcETransportTrap(void);
	mcETransportTrap(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);

	void beginTransport(mcParticle& p) override;
	void beginTransportInside(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
};

