// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <vector>

/// <summary>
/// Класс для поддержки самплинга параметров из комулятивной гистограммы в интервале аргументов от 0 до 1.
/// </summary>
class mcHistogramSampler
{
public:
	mcHistogramSampler();
	~mcHistogramSampler();

	/// <summary>
	/// Иницилизация данных из произвольного распределения, представленного бинами постоянного шага.
	/// </summary>
	/// <param name="np">количество точек в создаваемом самплере</param>
	/// <param name="nbins">количество бинов в исходном распределении</param>
	/// <param name="bin_size">начальная координата нулевого бина</param>
	/// <param name="bin_size">толщина (шаг) бина</param>
	/// <param name="data">массив распределения - количество самплов бинах</param>
	void setFromDistribution(int np, int nbins, double par_min, double bin_size, const double* data);
	void setFromDistribution(int np, double par_min, double bin_size, const std::vector<double>& data);

	/// <summary>
	/// Восстановление данных из буфера памяти.
	/// </summary>
	/// <param name="np">количество точек в самплере</param>
	/// <param name="data">массив данных готовой к употреблению гистограммы</param>
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
