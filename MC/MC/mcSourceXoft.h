// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// Класс сферического источника C60 
class mcSourceXoft : public mcSource
{
public:
	mcSourceXoft(const char* name, int nThreads,
		const geomVector3D& p, const geomVector3D& v, double r);
	virtual ~mcSourceXoft(void);

	void sample(mcParticle& p, mcThread* thread) override;

	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceXoft& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \tC60" << endl;
		os << "NAME = \t" << s.getName() << endl;
		os << "DIAMETER = \t" << 2 * s.r_ << endl;
		os << "POSITION = \t" << s.p_ << endl;
		return os;
	}

protected:
	geomVector3D p_;
	geomVector3D v_;
	double r_;
};