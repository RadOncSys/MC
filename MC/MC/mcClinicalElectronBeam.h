// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// Источник, представляющий реалистичный пучок облучения электронами.
// Указываемая позиция относится к изоцентру аппарата.
// Виртуальный источник предполагается находящимся в положении -100 см от него.
// Рамка поля находится в положении -5 см.
// В параметрах конструктора указыаются угловой разброс (сигма гауссиана) в радианах и энергетический разброс в МэВ.

class mcClinicalElectronBeam : public mcSource
{
public:
	mcClinicalElectronBeam(const char* name, int nThreads, double ke, double z0, 
		double sigmaE, double sigmaA, double fsx, double fsy);

	void sample(mcParticle& p, mcThread* thread) override;

	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcClinicalElectronBeam& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.p_ << endl;
		os << "DIRECTION = \t" << s.v_ << endl;
		os << "SIGMA_E = \t" << s.sigmaE_ << "\tradian" << endl;
		os << "SIGMA_A = \t" << s.sigmaA_ << "\tMeV" << endl;
		os << "FSX = \t" << s.fsx_ << endl;
		os << "FSY = \t" << s.fsy_ << endl;
		return os;
	}

protected:
	mc_particle_t type_;
	double ke_;
	geomVector3D p_;
	geomVector3D v_;
	int q_; // заряд
	double sigmaE_;
	double sigmaA_;
	double fsx_;
	double fsy_;
};
