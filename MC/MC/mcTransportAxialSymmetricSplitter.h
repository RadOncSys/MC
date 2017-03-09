// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

//  ласс транспорта типа специального фильтра.
// ¬ объекте этого класса не происходит взиамодействи€.
// ѕросто частицы, пересекающие плоскость в положительном направлении расщепл€ютс€
// на множество копий с пропорционально уменьшенными весами.
//  ажда€ копи€ поворачиваетс€ на свой угол вокруг оси Z
class mcTransportAxialSymmetricSplitter : public mcTransport
{
public:
	mcTransportAxialSymmetricSplitter(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, mc_particle_t ptype, int nsplit);
	~mcTransportAxialSymmetricSplitter(void);

	// Ќачало транспорта переписано, так как не нужно запускать симул€цию частиц.
	// »х просто нужно зарегистрировать и пропустить дальше.
	void beginTransport(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;

protected:
	enum mc_particle_t ptype_;
	int nsplit_;
};

