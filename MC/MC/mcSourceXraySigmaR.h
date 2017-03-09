// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// »сточник ускорител€ с вращающимс€ пучком и имеющим конечный размер,
// описываемы параметром сигма гаусового распределени€ электронов по радиусу.
// ¬ращение означает определенный угол падени€ ускоренных частиц на мишень.
class mcSourceXraySigmaR : public mcSource
{
public:
	mcSourceXraySigmaR(const char* name, int nThreads, mc_particle_t type, double ke, double z, double r, double theta);
	virtual ~mcSourceXraySigmaR(void);

	void sample(mcParticle& p, mcThread* thread) override;
	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceXraySigmaR& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.z_ << endl;
		os << "ROTATION ANGLE = \t" << s.theta_ << endl;
		return os;
	}

protected:
	mc_particle_t type_;
	int q_;			// зар€д
	double ke_;		// энерги€
	double z_;		// положение плоскости задани€ чатиц
	double sigma_;	// сигма разброса по одной оси
	double theta_;	// угол падени€ электронов на мишень (в градусах), обусловленный вращением

	// ¬спомогательные переменные
	double tr_;		// угол падени€ в радианах
	double uz_;		// проекци€ вектора скорости на ось z
	double sinu_;	// синус угла наклона вектора электронов
};
