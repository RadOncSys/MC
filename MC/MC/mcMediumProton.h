// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2021] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcMedium.h"
#include "mcDefs.h"

// Класс описания параметров конкретной среды для транспорта протонов.
class mcMediumProton : public mcMedium
{
public:
  mcMediumProton(void);
  virtual ~mcMediumProton(void);

  const double kEmax(void)const{return (double)dedx1_proto.size();}
  virtual void read(istream& is);
	const double AtomicWeight() const;	// Атомный вес среды, г/моль

private:
//--------------------------------
// Генерация данных (физика!)
//--------------------------------
	const double	gdEdxStragglingGaussVarianceConstPart(); // генерирует и возвращает постоянную (по энергии и пути) часть вариации Гаусссова приближения разброса dE/dx. 
	const double	gRadiationLength();	// генерирует и возвращает величину обратную радиационной длине сложного вещества rpp-2006-book.pdf 27.4.1 p.263 (eq.27.23) для расчёта радиационной длины отдельного элемента вызывает InverseRadiationLength  (принятое приближение (в версии 2007 года это приближение Dahl'а))
	const void		gSigmaInelastic(int Ap=1, int Zp=1);	// генерирует величину сечения неупругого взаимодействия

public:
  // TODO !!!
  // Надо разобраться с генерацией таблиц в идеологии PEGS
  // из таблиц тормозных способностей и физики протонов.
  // Для точного понимания смысла и форматов выходных данных
  // смотри класс mcPhysicsProton.

	//vector<double> dEdx_; // dE/dx - считываются из файла
	//vector<double> sigma_in_;	// сечение неупругих ядерных взаимодействий для энергий 1,2...,
	
	// не зависимая от энергии и пути часть Гауссовой вариации (sigma^2) dE/dx
	double dEdxStragglingGaussVarianceConstPart_;
	
	
  //// Обратить внимание, что в EGS в таблицах энергия задана в логарифмическом масштабе.
  //// Почему и нужно ли то же для протонов, если таблицы протонных dE/dX в линейном масштабе?

  //// Log energy index:
  //double iLogKE0_proto;
  //double iLogKE1_proto;

  //// Transport:
  double radLength;          // Radiation length, [cm](!)
  //double pStep;

  //// Moliere multiple scattering:
  //double scaledLittleB;
  //double chi_cc;

  //// Interactions:
  //// В случае электронов коэффициенты Branching определяли
  //// относительные вероятности двух моделируемых дискретных событий.
  //vector<double> br10_proto;         // Branching ratio 1
  //vector<double> br11_proto;

  //// Transport:
  // Коэффициенты интерполяции S(E)=S0[Ei]+S1[Ei]*E, Ei=int(E)
  // В электронах вместо энергии используется логарифм
  // Мы пока оставим линейную интерполяцию
  vector<double> sigma0_proto;       // Fictitious cross section, 1/cm
  vector<double> sigma1_proto;
  vector<double> dedx0_proto;        // Linear energy loss rate (stopping power)(Mev/cm)
  vector<double> dedx1_proto;
  //vector<double> stepSize0_proto;    // Step size
  //vector<double> stepSize1_proto;

  double transCutoff_proto; // Energy cutoff for proton transport
};
