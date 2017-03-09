// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <vector>

/// <summary>
/// ����� ��� ��������� ��������� ���������� �� ������������ ����������� � ��������� ���������� �� 0 �� 1.
/// </summary>
class mcHistogramSampler
{
public:
	mcHistogramSampler();
	~mcHistogramSampler();

	/// <summary>
	/// ������������ ������ �� ������������� �������������, ��������������� ������ ����������� ����.
	/// </summary>
	/// <param name="np">���������� ����� � ����������� ��������</param>
	/// <param name="nbins">���������� ����� � �������� �������������</param>
	/// <param name="bin_size">��������� ���������� �������� ����</param>
	/// <param name="bin_size">������� (���) ����</param>
	/// <param name="data">������ ������������� - ���������� ������� �����</param>
	void setFromDistribution(int np, int nbins, double par_min, double bin_size, const double* data);
	void setFromDistribution(int np, double par_min, double bin_size, const std::vector<double>& data);

	/// <summary>
	/// �������������� ������ �� ������ ������.
	/// </summary>
	/// <param name="np">���������� ����� � ��������</param>
	/// <param name="data">������ ������ ������� � ������������ �����������</param>
	void restore(int np, const double* data);
	void restore(int np, const std::vector<double>& data);

	unsigned npnt() const { return np_; }
	const double* data() const { return data_; }

	double sample(double f) const;

	friend std::ostream& operator << (std::ostream&, const mcHistogramSampler&);
	friend std::ostream& operator << (std::ostream&, const std::vector<mcHistogramSampler>&);

protected:
	void init(int np);
	void make_comulative(std::vector<double>& c, int nbins, const double* data);
	void make_comulative(std::vector<double>& c, const std::vector<double>& data);
	void set_from_comulative(int np, double par_min, double bin_size, const std::vector<double>& data);

protected:
	int np_;
	double step_;
	double* data_;
};
