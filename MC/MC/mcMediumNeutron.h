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

    const double kEmax(void)const{return (double)dedx1_proto.size();}
	const double AtomicWeight() const;

private:
	const double	gdEdxStragglingGaussVarianceConstPart(); 
	const double	gRadiationLength();
	const void		gSigmaInelastic(int Ap = 1, int Zp = 1);

public:
	double dEdxStragglingGaussVarianceConstPart_;
	double radLength;

	vector<double> sigma0_proto;
	vector<double> sigma1_proto;
	vector<double> dedx0_proto;
	vector<double> dedx1_proto;

	double transCutoff_neutron;
};
