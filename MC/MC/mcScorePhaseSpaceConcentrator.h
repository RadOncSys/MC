// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcscore.h"

/// <summary>
/// ����� ����� ���������� ������ ������ � ���� ���� ������������� � ��������� ������������ �������.
/// </summary>
class mcScorePhaseSpaceConcentrator : public mcScore
{
public:
	/// <summary>
	/// ����������� ������ ����� ����������� �� ������������� ������ �� ������� ��������� � �����
	/// � ��������������������� �������.
	/// </summary>
	/// <param name="module_name">��� ������, � ������� ������������ ������ scoring</param>
	/// <param name="nThreads">���������� �������</param>
	/// <param name="ptype">��� ������, ��� ������� ���������� ����������</param>
	/// <param name="focus">���������� �� ������������ ������, �� �������� ���������������� ������� ������� (��������, ������ ���������)</param>
	/// <param name="isXray">����� �������� ������, ��� �����-�������� � �������������� ��������� ������������ �������</param>
	/// <param name="ne">���������� ����� �� �������</param>
	/// <param name="nr">���������� ����� �� �������</param>
	/// <param name="nthet">���������� ����� �� ����</param>
	/// <param name="naxial">���������� ����� �� �������</param>
	/// <param name="emax">������������ ������� � ���</param>
	/// <param name="rmax">������������ ������ � ��</param>
	/// <param name="thetmax">������������ ���� � ��������</param>
	/// <param name="stat_file">��� ����� ������ �����������</param>
	mcScorePhaseSpaceConcentrator(
		const char* module_name, int nThreads, 
		enum mc_particle_t ptype, bool isXray, double focus,
		int ne, int nr, int nthet, int naxial, 
		double emax, double rmax, double thetmax, 
		const char* stat_file);
	virtual ~mcScorePhaseSpaceConcentrator();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScorePhaseSpaceConcentrator&);

protected:
	// ���������� ������� 4D ������� ����������, ���� ������ ���� ���������������� �������.
	int particleIdx(const mcParticle& particle) const;

	// ����� ���������� �� ���� �������������� sude by side
	// (��������������� � ������ ������� ��� ��������� ��������������, � ��������������� �������)
	void dumpDistributions(ostream& os, double* dsim, double* dmodel) const;

protected:
	enum mc_particle_t ptype_;
	int isxray_;
	int ne_;
	int nr_;
	int nthet_;
	int naxial_;

	double focus_;
	double emax_;
	double rmax_;
	double thetmax_;		// (1 - cos(thetmax_abs_ * Pi / 180))
	double thetmax_abs_;	// ������������ ���� � ��������

	std::string model_file_;

	// ��������� ����������
	int nne_;
	int nnr_;
	int nnt_;
	double de_;
	double dr_;
	double dthet_;
	double daxial_;
	double* data_;	// 4-� ������ ������ ������ ���������������� �����������
};
