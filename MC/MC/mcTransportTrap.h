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
class mcTransportTrap : public mcTransport
{
public:
	mcTransportTrap(void);
	mcTransportTrap(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportTrap(void);

	virtual void beginTransport(mcParticle& p);

	virtual void dumpVRML(ostream& os)const;

};
