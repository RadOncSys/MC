// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

class mcScoreMatrix2D : public mcScore
{
public:
	mcScoreMatrix2D(const char* module_name, int nThreads, int nx, int nz, const geomVector3D& v1, const geomVector3D& v2);
	virtual ~mcScoreMatrix2D(void);

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

	void CE2D() override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

public:
	int    Nx() const { return m_nx; }
	int    Nz() const { return m_nz; }
	double MinX() const { return m_xmin; }
	double MaxX() const { return m_xmax; }
	double MinY() const { return m_ymin; }
	double MaxY() const { return m_ymax; }
	double MinZ() const { return m_zmin; }
	double MaxZ() const { return m_zmax; }
	double StepX() const { return m_xstep; }
	double StepZ() const { return m_zstep; }
	void   mul(double f);
	//const double* GetMatrixConst() const { return m_M; }
	//double*       GetMatrix() { return m_M; }

	double	Dose(int iThread, int ir, int iz) const;
	double	Dose(int ir, int iz) const;
	double	Sigma(int ir, int iz) const;

protected:
	int     m_nx, m_nz;
	double  m_xmin, m_xmax;
	double  m_ymin, m_ymax;
	double  m_zmin, m_zmax;
	double  m_xstep, m_zstep;
	double* m_MAll;
	double** m_M;
};
