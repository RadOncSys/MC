// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcEndfP.h"
#include "mcMedium.h"
#include "mcDefs.h"

// Класс описания параметров конкретной среды для транспорта протонов.
class mcMediumProton : public mcMedium
{
public:
	mcMediumProton(void);
	virtual ~mcMediumProton(void);

	double kEmax(void)const { return (double)dedx1_proto.size(); }
	virtual void read(istream& is);
	void createDB();
	const double AtomicWeight() const;	// Атомный вес среды, г/моль

private:
	//--------------------------------
	// Генерация данных (физика!)
	//--------------------------------
	double gdEdxStragglingGaussVarianceConstPart();	// генерирует и возвращает постоянную (по энергии и пути) часть вариации Гаусссова приближения разброса dE/dx. 
	double gRadiationLength();	// генерирует и возвращает величину обратную радиационной длине сложного вещества rpp-2006-book.pdf 27.4.1 p.263 (eq.27.23) для расчёта радиационной длины отдельного элемента вызывает InverseRadiationLength  (принятое приближение (в версии 2007 года это приближение Dahl'а))
	void gSigmaInelastic(int Ap = 1, int Zp = 1);	// генерирует величину сечения неупругого взаимодействия

public:
	// не зависимая от энергии и пути часть Гауссовой вариации (sigma^2) dE/dx
	double dEdxStragglingGaussVarianceConstPart_;
	
	//// Transport:
	double radLength;          // Radiation length, [cm](!)

	// Transport:
	// Коэффициенты интерполяции S(E)=S0[Ei]+S1[Ei]*E, Ei=int(E)
	// В электронах вместо энергии используется логарифм
	// Мы пока оставим линейную интерполяцию
	vector<double> sigma0_proto;	// Fictitious cross section, 1/cm (речь о ядерных взаимодействиях)
	vector<double> sigma1_proto;
	vector<double> dedx0_proto;		// Linear energy loss rate (stopping power)(Mev/cm)
	vector<double> dedx1_proto;

	double transCutoff_proto;		// Energy cutoff for proton transport

	std::shared_ptr<std::vector<std::shared_ptr<mcEndfP>>> ENDFdata;
	double microsigmaforelement(int A, int Z, double kE) const;
};
