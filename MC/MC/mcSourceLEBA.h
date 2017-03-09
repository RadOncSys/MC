// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// Типы формы энергетическог спектра
enum spectrum_distr_t { SPECTRUM_GAUSS = 0, SPECTRUM_TRIANGLE, SPECTRUM_PRISM };

// Типы радиальных профилей
enum profile_distr_t { PROFILE_PRISM = 0, PROFILE_GAUSS, PROFILE_EXPONENT };

// Источник ускорителя со свойствами по мотивам симуляции Сергея ползова.
// От отличается экспоненциальным распределением интенсивности по радиусу
// и широким энергетическим спектром в форме прямоугольника или комбинации из двух треугольников.
// В данном источнике возможет только один треугольник.
// Угловой разброс частиц игнорируется из за малой величины.
// LEBA - Linear Electron Beam Accelerator.
class mcSourceLEBA : public mcSource
{
public:
	// name - нзвание модуля;
	// sptype - тип формы спектра (0 - гауссов, 1 - треугольный, 2 - прямоугльный);
	// kemean - средняя энергия в МэВ (середина прямоугольного спектра или пик других форм);
	// kewidth -ширина энергетического спектра в МэВ (ширина прямоугольника, ширина основания треугольника, или 2 сигма гауссиана);
	// z - положение истчника в см;
	// rsigma - показатель экспоненты распределения интенсивности по радиусу.
	mcSourceLEBA(const char* name, int nThreads, spectrum_distr_t sptype, profile_distr_t prf_type,
		double kemean, double kewidth, double z, double rsigma);

	void sample(mcParticle& p, mcThread* thread) override;
	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceLEBA& s)
	{
		os << (const mcSource&)s;
		os << "POSITION = \t" << s.z_ << endl;
		os << "TYPE = \t" << s.sptype_ << endl;
		os << "KEMEAN = \t" << s.kemean_ << endl;
		os << "KEWIDTH = \t" << s.kewidth_ << endl;
		os << "RSIGMA = \t" << s.rsigma_ << endl;
		return os;
	}

protected:
	spectrum_distr_t sptype_;
	profile_distr_t prf_type_;
	int q_;				// заряд
	double kemean_;		// средняя энергия
	double kewidth_;	// ширина энергетического спектра
	double z_;			// положение плоскости задания чатиц
	double rsigma_;		// показатель экспоненты распределения по радиусу
};
