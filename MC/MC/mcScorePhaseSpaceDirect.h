// Radiation Oncology Monte Carlo open source project
//
// Author: [2020] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

/// <summary>
/// ����� ����� ������ ������ � ��������� ������������ �������.
/// ����������� ��������� ������������ ������ ���������, ��������� �� ������ ������, 
/// ����������� ��������� ������ ��������� ����� ��������������������� ��������.
/// </summary>
class mcScorePhaseSpaceDirect : public mcScore
{
public:
	/// <summary>
	/// ����������� ������ ����� ����������� �� ������������� ������ �� ������� ��������� � �����
	/// � ��������������������� �������.
	/// </summary>
	/// <param name="module_name">��� ������, � ������� ������������ ������ scoring</param>
	/// <param name="nThreads">���������� �������</param>
	/// <param name="ptype">��� ������, ��� ������� ���������� ����������</param>
	/// <param name="emax">������������ ������� � ���</param>
	/// <param name="rmax">������������ ������ � ��</param>
	/// <param name="stat_file">��� ����� ������ �����������</param>
	mcScorePhaseSpaceDirect(const char* module_name, int nThreads, enum mc_particle_t ptype, 
		double emax, double rmax, const char* stat_file);
	virtual ~mcScorePhaseSpaceDirect();

	void ScoreFluence(const mcParticle& particle) override;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScorePhaseSpaceDirect&);

protected:
	enum mc_particle_t ptype_;
	double emax_;
	double rmax_;
	std::string model_file_;

	double energyScale_;
	double rScale_;				// �������� �������������� ��������
	double radialAngleScale_;	// ���������� Ux �������� �������� ([-1,1] -> [0,2])
	double azimutAngleScale_;	// ���������� Uy �������� �������� ([-1,1] -> [0,2])

	// ������� ���������� ������ (��� ��������� - threads, particles, parameters)
	unsigned short* data_;
	std::vector<unsigned> threadIndexes_;

	// Debug
	unsigned nr_;
	double dr_;
	std::vector<double> weight_;
	std::vector<double> fluence_;
	std::vector<double> rangle_;
	std::vector<double> aangle_;
	std::vector<double> rangle2_;
	std::vector<double> aangle2_;
};
