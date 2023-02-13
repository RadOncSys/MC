// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2022] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcMedium.h"
#include "mcDefs.h"

// Класс описания параметров конкретной среды для транспорта нейтронов.
class mcMediumNeutron : public mcMedium
{
public:
	mcMediumNeutron(void);
	virtual ~mcMediumNeutron(void);

	virtual void read(istream& is);

    const double kEmax(void)const{return (double)dedx1_neutro.size();}
	const double AtomicWeight() const;

private:
	const double	gdEdxStragglingGaussVarianceConstPart(); 
	const double	gRadiationLength();
	const void		gSigmaInelastic(int Ap = 1, int Zp = 1);

public:
	double dEdxStragglingGaussVarianceConstPart_;
	double radLength;

	vector<double> sigma0_neutro;
	vector<double> sigma1_neutro;
	vector<double> dedx0_neutro;
	vector<double> dedx1_neutro;

	double transCutoff_neutron;
};
