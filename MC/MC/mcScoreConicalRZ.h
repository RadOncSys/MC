// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

/// <summary>
/// Скоринг дозовых распределений в геометрии RZ с веерными гранями вдоль оси пучка.
/// </summary>
class mcScoreConicalRZ : public mcScore
{
public:
	/// <summary>
	/// Основной конструктора скоринга в веерной RZ геометрии
	/// </summary>
	/// <param name="module_name">произвольное название объекта</param>
	/// <param name="nThreads">количество обслуживаемых потоков</param>
	/// <param name="nr">количество колец скоринга</param>
	/// <param name="nz">количество слоев скоринга</param>
	/// <param name="rmax">радиус внешнего кольца на уровне изоцентра</param>
	/// <param name="zmin">глубина залегания передней поверхности скоринга в родительском фантоме</param>
	/// <param name="zmax">глубина залегания дальней поверхнсти скоринга</param>
	/// <param name="ziso">глубина залегания изоцентра в родительском фантоме</param>
	/// <param name="sad">расстояние - источник - изоцентр аппарата</param>
	mcScoreConicalRZ(const char* module_name, int nThreads, int nr, int nz, double rmax,
		double zmin, double zmax, double ziso, double sad);

	virtual ~mcScoreConicalRZ();

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

	friend ostream& operator << (ostream&, const mcScoreConicalRZ&);

public:
	int	Nr() const { return m_nr; }
	int	Nz() const { return m_nz; }
	double	MaxR() const { return m_rmax; }
	double	MinZ() const { return m_zmin; }
	double	MaxZ() const { return m_zmax; }
	double	StepR() const { return m_rstep; }
	double	StepZ() const { return m_zstep; }

	double	Dose(int iThread, int ir, int iz) const;
	double	Dose(int ir, int iz) const;
	double	Sigma(int ir, int iz) const;
	//double Dosemean (int iThread, int ir, int iz) const; средняя доза по вокселям
	//double Sigma (int iThread, int ir, int iz) const;

protected:
	int     m_nr, m_nz;     // число вокселов по радиусу (m_nr) и по z (m_nz)
	double  m_rmax, m_rm2;  // m_rm2 - квадрат максимального радиуса
	double  m_zmin, m_zmax;
	double  m_rstep, m_zstep;
	double  m_iso;
	double  m_sad;
	double* m_MAll;
	double** m_M;
};
