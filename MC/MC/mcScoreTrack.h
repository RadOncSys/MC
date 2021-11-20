// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcscore.h"
#include <vector>

// ����� ����������� ���������� ������ ����� �������� ����������.
class mcScoreTrack : public mcScore
{
public:
	mcScoreTrack(int nThreads, double R, double Z1, double Z2, double EMIN, bool doPhotons, bool doElectrons, bool doPositrons, bool doProtons, bool doNeutrons);
	~mcScoreTrack(void);

	void score(int iThread
		, mc_particle_t pt
		, const geomVector3D& p0
		, const geomVector3D& p1
		, double ke);

	void dumpVRML(ostream&) const override;

protected:
	vector<vector<geomVector3D>> photons_;
	vector<vector<geomVector3D>> electrons_;
	vector<vector<geomVector3D>> positrons_;
	vector<vector<geomVector3D>> protons_;
	vector<vector<geomVector3D>> neutrons_;

	// ������� ��������, �������������� ����������� ������
	double R_, Z1_, Z2_;

	// ����������� �������, ���� ������� ����� �� ������������
	double EMIN_;

	// ����� ������������ ������
	bool doPhotons_, doElectrons_, doPositrons_, doProtons_, doNeutrons_;
};
