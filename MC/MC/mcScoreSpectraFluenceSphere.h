// Radiation Oncology Monte Carlo open source project
//
// Author: [2020] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// ����� �������� �������� ������ ��������� � ����������� �������.
// ���������� ������� �� ������� � ������� ����������� ������� �� �� ���������� (����������� �� ���).
// ����������� �������������� ������� ���� �� �������, �� ������������� ����� �������������� �������
// (����� ��������� ��������� ������������� �� ������).
// ����������� ���������� ����� ������� � ��� �� ���������� �����.
// ����������� ������ ������� ���������� ���� � ������� ������������ ������� ���� ���������� ������.
// ��� ������ ���������� �������� ������� �� ������� ������ � ��^2.
class mcScoreSpectraFluenceSphere : public mcScore
{
public:
	mcScoreSpectraFluenceSphere(const char* module_name, int nThreads, mc_particle_t pt, 
		double ecut, int nr, int ne, double rmax, double emax, double radius);
	virtual ~mcScoreSpectraFluenceSphere();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreSpectraFluenceSphere&);

protected:
	double eFluenceBin(int ir) const;
	double nFluenceBin(int ir) const;
	double nSpectrumBin(int ir, int ie) const;
	double eSpectrumBin(int ir, int ie) const;
	double sAngle(int ir) const;

	mc_particle_t pt_;
	int nr_;
	double rstep_;	// ��� �� ������� � ��������
	double rmax_;
	int ne_;
	double estep_;
	double emax_;
	double ecut_;	// �������, ���� ������� ������� ����������� �� ����������
	double radius_;	// Detector radius (for VRML only)
	std::vector<std::vector<double>> energy_fluence_;
	std::vector<std::vector<double>> number_fluence_;
	std::vector<std::vector<std::vector<double>>> number_spectra_;
	std::vector<std::vector<std::vector<double>>> energy_spectra_;
};
