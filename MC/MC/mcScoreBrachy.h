// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

// —коринг брахитерапевтического источника.
// Ќаиболее близок к mcScoreMatrixRZ.
// ¬ отлиичие от последнего не выводит дозовую матрицу из-за ее большого размера.
// ¬место этого в статистике выводит сразу даблицы g(r) и F(r,thetta).
//  роме того, размеры матрицы не читаютс€ их конфигурации расчетов, а устанавливаютс€ стстически.
// ћатрица F(r,thetta) благодар€ удалению обратных кадратов мен€етс€ не резко.
// ѕоэтому она после расчета подвергаетс€ сглаживанию.

class mcScoreBrachy : public mcScore
{
public:
	mcScoreBrachy(const char* module_name, int nThreads, double sourceLength);

	virtual ~mcScoreBrachy();

	void ScoreFluence(const mcParticle& particle) override;

	void ScorePoint(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0) override;

	void ScoreLine(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0
		, const geomVector3D& p1) override;

	virtual void scoreEnergyInVoxel(int iThread, int ir, int iz, double edep);

	void CE2D() override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreBrachy&);

public:
	//int	Nr() const { return m_nr; }
	//int	Nz() const { return m_nz; }
	//double	MaxR() const { return m_rmax; }
	//double	MinZ() const { return m_zmin; }
	//double	MaxZ() const { return m_zmax; }
	//double	StepR() const { return m_rstep; }
	//double	StepZ() const { return m_zstep; }

	double	Dose(int iThread, int ir, int iz) const;
	double	Dose(int ir, int iz) const;

protected:
	int m_nr, m_nz;			// число вокселов по радиусу (m_nr) и по z (m_nz)
	double m_rmax, m_rm2;	// m_rm2 - квадрат максимального радиуса
	double m_zmin, m_zmax;
	double m_rstep, m_zstep;
	double L_;
	double* m_MAll;
	double** m_M;
};
