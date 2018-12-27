// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// ����� ��������� �������� ������������� � ����������� ������.
// ���������� ������������ ��� ������� ������������ ������ ������� �������������������� ���������.
// �������� ����� ������ ���� ���� �� �������, �� �� ����������� ��� ��������� �������
// (�.�. ������������� �������������� ���������������� ��-�� ������������� ������.
// ������� �� ������ ��� �� ������� �� ������� � ����������.
// �������� ���������� �������� � ������� (��� Z) ��������� �� ������ ���������� 
// ������� ����������� ������� � ���������� ������� �� ���������.
class mcScoreSphereMatrix : public mcScore
{
public:
	mcScoreSphereMatrix(const char* module_name, int nThreads, int np, int nm, double rmin, double rmax);
	virtual ~mcScoreSphereMatrix();

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

	virtual void scoreEnergyInVoxel(int iThread, int ip, int im, double edep);

	void CE2D() override;

	double	Dose(int ip, int im) const;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreSphereMatrix&);

protected:
	// ���������� �� ����� p0 � ����������� v �� ������, ���������������� ������ angle
	double distanceToParallel(const geomVector3D& p0, const geomVector3D& v, int ip, bool isPOutside);

	// ���������� �� ����� p0 � ����������� v �� ���������, ��������������� ������� angle
	double distanceToMeridian(const geomVector3D& p0, const geomVector3D& v, int im, bool isMPositive);

protected:
	int np_;				// ���������� ����������
	int nm_;				// ���������� ����������
	double rmin_, rmax_;	// ������� ����, ����� �������� ���� ���������
	double* MAll_;
	double** M_;

	double r2min_, r2max_;
	double pstep_;
	double mstep_;
};
