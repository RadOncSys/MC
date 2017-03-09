// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта типа специального фильтра.
// В объекте этого класса не происходит взиамодействие.
// Просто частицы, пересекающие плоскость в положительном направлении расщепляются
// на множество копий с пропорционально уменьшенными весами.
class mcTransportSimpleSplitter : public mcTransport
{
public:
	mcTransportSimpleSplitter(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, mc_particle_t ptype, int nsplit);
	~mcTransportSimpleSplitter(void);

	// Начало транспорта переписано, так как не нужно запускать симуляцию частиц.
	// Их просто нужно зарегистрировать и пропустить дальше.
	void beginTransport(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;

protected:
	enum mc_particle_t ptype_;
	int nsplit_;
};
