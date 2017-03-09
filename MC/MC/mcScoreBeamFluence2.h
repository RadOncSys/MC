// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

class mcScoreBeamFluence2 : public mcScore
{
public:
	mcScoreBeamFluence2(const char* module_name, int nThreads, int nr, double rmax, int nr_s, double rmax_s, double H, int nz, double d_iz);
	virtual ~mcScoreBeamFluence2();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreBeamFluence2&);

public:
	double Intencity(int iThread, int ir) const;
	double Intencity(int ir) const;
	double MeanEnergy(int iThread, int ir) const;
	double MeanEnergy(int ir) const;

	double Intencity_s(int iThread, int iz, int ir) const;
	double Intencity_s(int iz, int ir) const;
	double Value_s(int iz, int ir) const;

protected:
	int     m_nr;
	double  m_rmax;
	double  m_rstep;
	double* m_intencity_all;
	double** m_intencity;
	double* m_intencity_count_all;
	double** m_intencity_count;

	double	h_;	// растояние от скоринга до плоскости источника
	int		nz_; // количесво слоёв скоринга по оси oZ
	double	d_iz_; // расстояние между слоями скоринга
	int     m_nr_s;
	double  m_rmax_s;
	double  m_rstep_s;
	double*** m_intencity_s;
};
