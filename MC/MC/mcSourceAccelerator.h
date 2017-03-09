// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// Источник ускорителя с вращающемся пучком и
// имеющим конечный размер, описываемы параметром радиус.
// Вращение означает определенный угол падения ускоренных частиц на мишень.
// Предполагается, что поток частиц равнмерно рспределен в сечении пучка.
class mcSourceAccelerator : public mcSource
{
public:
	mcSourceAccelerator(const char* name, int nThreads, mc_particle_t type, double ke, double z, double r, double theta);

	void sample(mcParticle& p, mcThread* thread) override;
	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceAccelerator& s)
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
	int q_;			// заряд
	double ke_;		// энергия
	double z_;		// положение плоскости задания чатиц
	double r_;		// радиус пучка
	double theta_;	// угол падения электронов на мишень (в градусах), обусловленный вращением

	// Вспомогательные переменные
	double tr_;		// угол падения в радианах
	double uz_;		// проекция вектора скорости на ось z
	double sinu_;	// синус угла наклона вектора электронов
};
