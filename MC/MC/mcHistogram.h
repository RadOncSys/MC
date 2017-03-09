// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <vector>

/// <summary>
/// Структура, обслуживающая гистограммы. В экземпляре содержится точка для добавления в гистограмму.
/// </summary>
struct mcHPair
{
	mcHPair() :d_(0), w_(0) {};
	mcHPair(double d, double w) :d_(d), w_(w) {};
	mcHPair(double d) :d_(d), w_(1.) {};
	double d_;
	double w_;
};

/// <summary>
/// Класс гистограмм для поддержки .
/// Данный источник дергает частицы разыгрывая параметры из интегральных диаграм,
/// полученных из реального фазового пространства.
/// </summary>
class mcHistogram
{
public:
	mcHistogram();
	mcHistogram(double min, double step, int np);
	mcHistogram(const mcHistogram&);

	unsigned npnt() const { return np_; }
	unsigned nmax() const { return nmax_; }
	double argmin() const { return min_; }
	double argmax() const { return max_; }
	double step() const { return step_; }
	double total() const { return total_; }
	double sum() const { return sum_; }
	double sum2() const { return sum2_; }
	double arg(unsigned i) const { return min_ + (double(i) + 0.5)*step_; }
	const std::vector<double>& d() const { return d_; }

	double mean()const;
	double median()const;
	double stdev()const;
	double sample(double f)const;

	void init(double min, double step, unsigned np);
	void normalize();
	void makePerArea();
	void makePerSolidAngle();
	void makeComulative();
	void operator+=(const mcHPair&);
	void operator+=(const mcHistogram&);
	friend std::ostream& operator << (std::ostream&, const mcHistogram&);
	friend std::ostream& operator << (std::ostream&, const std::vector<mcHistogram>&);

protected:
	bool comulative_;
	unsigned np_;
	unsigned nmax_;
	double min_;
	double max_;
	double step_;
	double total_;
	double sum_;
	double sum2_;
	std::vector<double> d_;
};
