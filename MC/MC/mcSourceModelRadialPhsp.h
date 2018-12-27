// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"
#include "mcHistogramSampler.h"

// ���� ������������, ����� �� ������������� �� ������� ����������������� ��������� 
// ������� ������� � ���������� ���������. 
// � ��������� ������ ��� ����� ����� ����������������� �������������.
// � �������� ����������� �������� ���������, ��� ��� ������������ ���������� ����� 
// ��� ������������� �� ����������� ������ ������. 
// �������, �� ������� ���� ��� ����� ���������� ����� ������������ ������� ������������ ����������.

//#define SRC_USE_ASIMUT_MODEL

/// <summary>
/// ����� ������ ��������� C60 �������������� ��������.
/// ������ �������� ������� ������� ���������� ��������� �� ������������ �������,
/// ���������� �� ��������� �������� ������������.
/// </summary>
class mcSourceModelRadialPhsp : public mcSource
{
public:
	mcSourceModelRadialPhsp(const char* name, int nThreads, double z0);

	void sample(mcParticle& p, mcThread* thread) override;
	void sample(mcParticle& p, mcRng& rng) const;

	void dumpVRML(ostream& os) const override;

	friend std::ostream& operator << (std::ostream&,const mcSourceModelRadialPhsp&);

	/// <summary>
	/// �������� ������ ������ � ������. ������ ������������� ������ �������.
	/// ������� ���������� ������� ������� ���������� �� � ������� free();
	/// �������� size �������� ����� ���������� ������ � ������
	/// </summary>
	void* saveToMemory(int& size);

	/// <summary>
	/// �������������� ������ �� ������ ������.
	/// </summary>
	void readFromMemory(void* buffer);

	/// <summary>
	/// �������� ��������� �� ������� ���������.
	/// </summary>
	void createFromDistribution(int isxray, int ne, int nr, int nt, int na,
		double emax, double r2max, double tmax, double focus,
		const double* data);

	static double sampleAzimut(const vector<double>& a, double x);
	static void fitAsimutDistribution(int np, const double* data, double& a, double& b);

protected:
	/// <summary>
	/// ��������� ��������� � ����� ��� �������� ������� ���������������� ����.
	/// </summary>
	void sampleForEnergy(int eidx, double& r, double& t, double& a, double rnd0, double rnd1, double rnd2) const;

	/// <summary>
	/// �������� �������� ���������.
	/// </summary>
	void createHistograms(int ne, int nr, int nt, int na);

protected:
	int isxray_;

	// ��������� ��������� ��������� �� ��� �������
	double z0_;

	int ne_;
	int nr_;
	int nthet_;
	int naxial_;

	double focus_;
	double emax_;
	double rmax_;
	double thetmax_abs_;
	double thetmax_;

	// ��������� �������������� ����������
	double de_;
	double dr_;
	double dthet_abs_;
	double dthet_;
	double daxial_;

	// ���� ����� �������
	double we117_;
	double we133_;
	double wexscale_;

	// ����������� ��� ��������� ��������� � �����
	mcHistogramSampler xehistogram_;
	vector<mcHistogramSampler> rhistograms_;
	vector<vector<mcHistogramSampler>> thistograms_;

#ifdef SRC_USE_ASIMUT_MODEL
	vector<vector<vector<vector<double>>>> apars_;
#else
	vector<vector<vector<mcHistogramSampler>>> ahistograms_;
#endif
};
