#include "mcPhysicsCommon.h"
#include "mcDefs.h"
#include "mcRng.h"
#include <math.h>

//=================================
// Standard Physical functions
//---------------------------------

// Возвращает энергию центра масс для покоящейся частицы c targetMass и 
// налетающей с projectileMass и ПОЛНОЙ энергией totalProjectileEnergy
// Всё в MeV
// rpp p.321, eq.38.3
double ECenterOfMass(double targetMass, double projectileMass, double totalProjectileEnergy) {
	double a = SQUARE(targetMass) + SQUARE(projectileMass) + 2 * targetMass*totalProjectileEnergy;
	a = sqrt(a);
	return a;
}

//=================================
// Основные релятивистские формулы
//---------------------------------

// Релятивистский фактор gamma=sqrt(1/(1-beta^2))
// tEnergy - полная энергия частицы, Mass - её масса
inline double gamma(double Mass, double tEnergy)
{
	return tEnergy / Mass;
}

// Скорость в единицах скорости света beta_squared=(velocity/light velocity)^2
// tEnergy - полная энергия частицы, Mass - её масса
double betasq(double Mass, double tEnergy) // 
{
	double	g_1 = Mass / tEnergy;
	return	1 - SQUARE(g_1);
}

// (beta*gamma) squared = P/M
// tEnergy - полная энергия частицы, Mass - её масса
double betagammasq(double Mass, double tEnergy)
{
	double	g = gamma(Mass, tEnergy);
	return	SQUARE(g) - 1;
}

//=================================
// Standard Random Numbers functions
//---------------------------------
// Гауссово распределённая случайная величина
// Получаем 2-е случайные величины r1 r2, стандартному нормальному закону
// (Гауссиана со средним = 0 и стандартным отклонением 1)
// методом Marsaglia
// Две потому что таково устройство метода
// см. http://en.wikipedia.org/wiki/Marsaglia_polar_method
// Генератор равномерно распределённых чисел rng.rnd() должен давать [0,1] -?
int GaussStandardRnd_by_Marsaglia(mcRng& rng, double& r1, double& r2)
{
	double x, y, s2;
	do {
		x = 2 * rng.rnd() - 1;
		y = 2 * rng.rnd() - 1;
		s2 = SQUARE(x) + SQUARE(y);
	} while (s2 >= 1);
	r1 = sqrt(-2 * log(s2) / s2);
	r2 = y*r1;
	r1 = x*r1;
	return 0;
}