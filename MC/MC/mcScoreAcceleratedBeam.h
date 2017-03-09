// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

/// <summary>
/// ����� ������� ����� ���������� ���������� �� ������ ������.
/// </summary>
class mcScoreAcceleratedBeam : public mcScore
{
public:
	/// <summary>
	/// ����������� ������ ������� ����� ���������� ���������� �� ������ ������.
	/// ���������� ��������� �� ���������������. ������, �� ������������ �������� ��������� �������.
	/// ������ ����� ���������� ������ ������� ��������, ��� ��� ������ �� ����� ���������
	/// � ��������� ���������� �������� ������������ ������.
	/// ���������� ���������� ������� �������� ����� ������ ��������� � ����� ����� � ����� � ���� � ������������ ��������.
	/// </summary>
	/// <param name="module_name">��� ������, � ������� ������������ ������ scoring</param>
	/// <param name="nThreads">���������� �������</param>
	/// <param name="ptype">��� ������, ��� ������� ���������� ����������</param>
	/// <param name="ne">���������� ����� �� �������</param>
	/// <param name="nr">���������� ����� �� �������</param>
	/// <param name="naxial">���������� ����� �� ������� ����������� �������</param>
	/// <param name="nthet">���������� ����� �� ����</param>
	/// <param name="emin">���������� ������� � ���</param>
	/// <param name="emax">������������ ������� � ���</param>
	/// <param name="rmax">������������ ������ � ��</param>
	/// <param name="thetmax">������������ ���� � ��������</param>
	mcScoreAcceleratedBeam(const char* module_name, int nThreads, enum mc_particle_t ptype,
		int ne, int nr, int naxial, int nthet,
		double emin, double emax, double rmax, double thetmax);

	void ScoreFluence(const mcParticle& particle) override;

	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreAcceleratedBeam&);

protected:
	enum mc_particle_t ptype_;
	int ne_;
	int nr_;
	int nthet_;
	int naxial_;

	double emin_;
	double emax_;
	double rmax_;
	double thetmax_;

	// ��������� ����������
	double de_;
	double dr_;
	double dthet_;
	double daxial_;

	double skipped_energy_;

	std::vector<double> espec_;				// ������ ��������� �� ���� ��������
	std::vector<double> rspec_;				// ������������� ������� �� ������� �� ���� ��������
	std::vector<double> aspec_;				// ������������� ������� �� ������� (������� ��������) �� ���� ��������
	std::vector<double> tspec_;				// ������������� ������� �� ����������� �� ���� ��������

	std::vector<vector<double>> evr_;		// �������������� ������� � ����������� �� �������
	std::vector<vector<double>> eva_;		// �������������� ������� � ����������� �� ������� (������� ��������)
	std::vector<vector<double>> tva_;		// ������� ������������� � ����������� �� �������
	std::vector<vector<double>> wvxy_;		// ����� ������� � ����������� �����������
	std::vector<vector<double>> wvvxvy_;	// ����� ������� � ����������� �� ����
};
