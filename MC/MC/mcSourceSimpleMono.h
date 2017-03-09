// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

class mcSourceSimpleMono : public mcSource
{
public:
	mcSourceSimpleMono(void);
	mcSourceSimpleMono(const char* name, int nThreads, mc_particle_t type, double ke, const geomVector3D& p, const geomVector3D& v);
	virtual ~mcSourceSimpleMono(void);

	void init(mc_particle_t type        // тип частиц
		, double ke                // кинетическая энергия
		, const geomVector3D& p    // точка рождения частиц
		, const geomVector3D& v);  // направление движения (единичный вектор)

	void sample(mcParticle& p, mcThread* thread) override;

	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceSimpleMono& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.p_ << endl;
		os << "DIRECTION = \t" << s.v_ << endl;
		return os;
	}

protected:
	mc_particle_t type_;
	double ke_;
	geomVector3D p_;
	geomVector3D v_;
	int q_; // заряд
};
