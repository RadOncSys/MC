
#pragma once
#include "mcsource.h"


class mcGaussSource : public mcSource
{
public:
	mcGaussSource(const char* name, int nThreads, mc_particle_t type, double ke,double spreadke, double z, double sigmax, double sigmay, double sigmathetax, double sigmathetay );
	virtual ~mcGaussSource(void);

	void sample(mcParticle& p, mcThread* thread) override;
	//void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcGaussSource& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.z_ << endl;
		/*os << "ROTATION ANGLE = \t" << s.theta_ << endl;*/
		return os;
	}

protected:
	mc_particle_t type_;
	int q_;			// заряд
	double ke_;		// энергия
	double z_;		// положение плоскости задания чатиц
	double sigmax_;	// сигма разброса по оси X
	double sigmay_; 	// сигма разброса по оси Y
	double sigmathetax_; // сигма разброса угла по оси X в радианах
	double sigmathetay_; // сигма разброса угла по оси Y в радианах
	double spreadke_; // равономерный разброс по энергии
	


	// Вспомогательные переменные
	//double tr_;		// угол падения в радианах
	//double uz_;		// проекция вектора скорости на ось z
	//double sinu_;	// синус угла наклона вектора электронов
};
