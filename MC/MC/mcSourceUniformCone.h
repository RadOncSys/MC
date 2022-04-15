// Radiation Oncology Monte Carlo open source project
//
// Author: [2022] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// Simple source of monoenergetic particles uniformly distributed within circular cone.
class mcSourceUniformCone : public mcSource
{
public:
	mcSourceUniformCone(const char* name, int nThreads, mc_particle_t type, double ke, double z, double r, double sad);
	virtual ~mcSourceUniformCone(void);

	void sample(mcParticle& p, mcThread* thread) override;
	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceUniformCone& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.z_ << endl;
		os << "RADIUS = \t" << s.r_ << endl;
		os << "SAD = \t" << s.sad_ << endl;
		return os;
	}

protected:
	mc_particle_t type_;
	int q_;
	double ke_;		
	double z_;		
	double r_;		// beam radius at distance sad_
	double sad_;	// sad_ is used as distance at which beam radius is defined and draw in wrml
};
