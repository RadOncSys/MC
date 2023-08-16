// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

//  ласс брахитерапевтического источника.
// Ёто просто цилиндр, внутри которого изотропно рождаютс€ частицы.
// ƒл€ симул€ции спектра используютс€ веса.
// Ќапример IR-192 представл€етс€ дес€тью энерги€ми фотонов, 
// по€вл€ющихс€ с одинаковой частотой но разными весам.
// ћатериал истчника здесь не важен, так как это проблема объекта транспорта
// ƒанный класс симулирует несколько типов источников (изотопов), указываемых типом
//---------------------------------------------------------------------------------

enum mc_isotope_t { C60 = 0, IR192, UNKNOWN };


class mcBrachySource : public mcSource
{
public:
	mcBrachySource(const char* name, int nThreads,
		const geomVector3D& p, const geomVector3D& v, double r, double h, mc_isotope_t isotope);
	virtual ~mcBrachySource(void);

	void sample(mcParticle& p, mcThread* thread) override;

	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcBrachySource& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \tC60" << endl;
		os << "NAME = \t" << s.getName() << endl;
		os << "DIAMETER = \t" << 2 * s.r_ << endl;
		os << "HEIGHT = \t" << s.h_ << endl;
		os << "POSITION = \t" << s.p_ << endl;
		return os;
	}

protected:
	mc_isotope_t isotope_;
	geomVector3D p_;
	geomVector3D v_;
	double r_;
	double h_;
};
