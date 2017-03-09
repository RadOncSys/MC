// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"

// ����� ������� ���������� � ��������� XY.
// ������� �����������.
// ��������� ������������ ����������� ������������ ��� X � Y.
/// Scoring �������� ������������� � ��������� XY
class mcScoreMatrixXY : public mcScore
{
public:
	/// <summary>
	/// ������������ �����������.
	/// </summary>
	/// <param name="module_name">��� ������ ����������, � �������� ������������� scoring.</param>
	/// <param name="nThreads">���������� ��������������� �������.</param>
	/// <param name="nx">������ ������� �� X.</param>
	/// <param name="ny">������ ������� �� Y.</param>
	/// <param name="psx">������ voxel X.</param>
	/// <param name="psy">������ voxel Y.</param>
	/// <param name="psz">������ voxel Z.</param>
	/// <param name="z0">������� ���������� �������� ����������� �������.</param>
	mcScoreMatrixXY(const char* module_name, int nThreads, int nx, int ny, double psx, double psy, double psz, double z0);
	virtual ~mcScoreMatrixXY(void);

	void ScoreFluence(const mcParticle& particle) override {}

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

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

public:
	double	Dose(int iThread, int ix, int iy) const;
	double	Dose(int ix, int iy) const;

protected:
	int    nx_, ny_;
	double z0_;
	double psx_;
	double psy_;
	double psz_;
	double* MAll_;
	double** M_;

	// ��������� ����������
	double z1_;
	double minx_, miny_;
};
