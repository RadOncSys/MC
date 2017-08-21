// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

class mcScoreMatrixCylinderAzimut : public mcScore
{
public:
	mcScoreMatrixCylinderAzimut(const char* module_name, int nThreads, int nr, int nz, double rmax, double zmin, double zmax);

	virtual ~mcScoreMatrixCylinderAzimut();

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

	virtual void scoreEnergyInVoxel(int iThread, int ir, int ia, double edep);

	void CE2D() override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreMatrixCylinderAzimut&);

public:
	int	Nr() const { return m_nr; }
	int	Na() const { return m_na; }
	double	MaxR() const { return m_rmax; }
	double	MinZ() const { return m_zmin; }
	double	MaxZ() const { return m_zmax; }
	double	StepR2() const { return m_r2step; }
	double	StepA() const { return m_astep; }

	double	Dose(int iThread, int ir, int iz) const;
	double	Dose(int ir, int iz) const;

protected:
	// ќпределение индексов €чейки дл€ указанной точки.
	// ќтрицательный индекс означает, что находимс€ за пределами скоинга.
	void getVoxelAtPoint(const geomVector3D& p, int& ir, int& ia);

	int     m_nr, m_na;     // число вокселов по радиусу (m_nr) и по азимуту (m_na)
	double  m_rmax, m_rm2;  // m_rm2 - квадрат максимального радиуса
	double  m_zmin, m_zmax;
	double  m_r2step, m_astep;
	double* m_MAll;
	double** m_M;
};
