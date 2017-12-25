// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

class mcHistogramSampler;

// Model of breamstrahlung radiation beam from classical electron accelerator.
// Simplified model takes uniform photon spectum (from commissioned TPS),
// samples uniformly distributed within the field rectangle,
// calculates particle direction and begin it's transport from the radiation source point.
// Note: It is legal to put lat say plexiglass slab at any distance from the radiation source 
// to mimic contamination electrons fluence from the radiation head.
class mcSourceXrayBeam : public mcSource
{
public:
	mcSourceXrayBeam(const char* name, int nThreads, double z, double enom, double sad, double x1, double x2, double y1, double y2);
	virtual ~mcSourceXrayBeam(void);

	void sample(mcParticle& p, mcThread* thread) override;
	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceXrayBeam& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "ENOM = \t" << s.nom_energy_ << endl;
		os << "SAD = \t" << s.sad_ << endl;
		os << "FSX1 = \t" << s.fsx1_ << endl;
		os << "FSX2 = \t" << s.fsx2_ << endl;
		os << "FSY1 = \t" << s.fsy1_ << endl;
		os << "FSY2 = \t" << s.fsy2_ << endl;
		return os;
	}

protected:
	mc_particle_t type_;
	double z_;			// radiation source position

	double nom_energy_;	// X-ry beam nominal energy
	double sad_;		// source - axis distance, at wich fiel size is defined
	double fsx1_;		// independent jaw X1 gjsition at sad
	double fsx2_;
	double fsy1_;
	double fsy2_;

	mcHistogramSampler* esampler_;	// energy sampler from spectrum
};
