// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcRng;

//=================================
// Standard Physical functions
//---------------------------------

// Возвращает энергию центра масс для покоящейся частицы c targetMass и 
// налетающей с projectileMass и ПОЛНОЙ энергией totalProjectileEnergy
// Всё в MeV
// rpp p.321, eq.38.3
double ECenterOfMass(double targetMass, double projectileMass, double totalProjectileEnergy);

//=================================
// Основные релятивистские формулы
//---------------------------------

// Релятивистский фактор gamma=sqrt(1/(1-beta^2))
// tEnergy - полная энергия частицы, Mass - её масса
inline double gamma(double Mass, double tEnergy);

// Скорость в единицах скорости света beta_squared=(velocity/light velocity)^2
// tEnergy - полная энергия частицы, Mass - её масса
double betasq(double Mass, double tEnergy);

// (beta*gamma) squared = P/M
// tEnergy - полная энергия частицы, Mass - её масса
double betagammasq(double Mass, double tEnergy);

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
int GaussStandardRnd_by_Marsaglia(mcRng& rng, double& r1, double& r2);
