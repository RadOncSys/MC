#include "mcHistogram.h"
#include <strstream>

using namespace std;

mcHistogram::mcHistogram()
	:comulative_(false)
	, np_(0)
	, nmax_(0)
	, min_(0)
	, max_(0)
	, step_(0)
	, total_(0)
	, sum_(0)
	, sum2_(0)
{
}

mcHistogram::mcHistogram(double min, double step, int np)
	:comulative_(false)
	, nmax_(0)
	, max_(0)
	, total_(0)
	, sum_(0)
	, sum2_(0)
	, d_(np, 0)
{
	init(min, step, np);
}

mcHistogram::mcHistogram(const mcHistogram& h)
	:comulative_(h.comulative_)
	, np_(h.npnt())
	, nmax_(h.nmax())
	, min_(h.argmin())
	, max_(h.argmax())
	, step_(h.step())
	, total_(h.total())
	, sum_(h.sum())
	, sum2_(h.sum2())
	, d_(h.d())
{
}

void mcHistogram::init(double min, double step, unsigned np)
{
	min_ = min;
	step_ = step;
	np_ = np;
	d_.resize(np_, 0);
}

void mcHistogram::normalize()
{
	if (comulative_)
		throw std::exception("попытка нормировки интегральной гистограммы");
	// Чтобы не было проблем из-за повторной нормировки сумму считаем заново
	unsigned i;
	double sum = 0;
	for (i = 0; i < np_; i++) sum += d_[i];
	if (sum == 0) return;
	for (i = 0; i < np_; i++) d_[i] /= sum;
}

double mcHistogram::mean()const
{
	return (total_ > 0) ? sum_ / total_ : 0.;
}

double mcHistogram::median()const
{
	if (comulative_)
		return sample(0.5);

	double f0 = 0, f1 = d_[0], f = 0.5;
	unsigned i;
	for (i = 1; i < np_; i++) {
		f1 += d_[i];
		if (f1 >= f) break;
		f0 = f1;
	}
	return (f1 == 0) ? 0 : (f0 == f1) ? arg(i - 1) : min_ + (i + 0.5 + (f - f0) / (f1 - f0))*step_;
}

double mcHistogram::stdev()const
{
	return (total_ > 0) ? sqrt((sum2_ - sum_*sum_ / total_) / total_) : 0.;
}

void mcHistogram::makePerArea()
{
	if (comulative_)
		throw std::exception("попытка нормировки интегральной гистограммы на площадь");
	double r1 = 0, r2;
	for (unsigned i = 0; i < np_; i++) {
		r2 = min_ + (i + 1)*step_;
		d_[i] /= 3.14*(r2*r2 - r1*r1);
		r1 = r2;
	}
}

void mcHistogram::makePerSolidAngle()
{
	if (comulative_)
		throw std::exception("попытка нормировки интегральной гистограммы на площадь");
	double r1 = 1, r2;
	for (unsigned i = 0; i < np_; i++) {
		r2 = cos(min_ + (i + 1)*step_);
		d_[i] /= 6.28*(r1 - r2);
		r1 = r2;
	}
}

void mcHistogram::makeComulative()
{
	if (comulative_)
		throw std::exception("повторная попытка преобразования в интегральную гистограмму");
	normalize();
	for (unsigned i = 1; i < np_; i++)
		d_[i] += d_[i - 1];
	comulative_ = true;
}

double mcHistogram::sample(double f)const
{
	if (!comulative_)
		throw std::exception("попытка самплинга дифференциальной гистограммы");
	if (f < 0 || f>1)
		throw std::exception("самплинг может осуществляться в интервале от 0 до 1");
	unsigned i;
	for (i = 1; i < np_ - 1; i++) if (d_[i] >= f) break;
	double d1 = d_[i - 1], d2 = d_[i];
	return (d1 == d2) ? arg(i - 1) : min_ + (i + (f - d2) / (d2 - d1))*step_;
}

void mcHistogram::operator+=(const mcHPair& pr)
{
	if (comulative_)
		throw std::exception("попытка добавления данных к интегральной гистограмме");
	if (pr.d_ < min_) return;
	unsigned i = unsigned((pr.d_ - min_) / step_);
	if (i >= np_) return;
	d_[i] += pr.w_;
	total_ += pr.w_;
	sum_ += pr.d_*pr.w_;
	sum2_ += pr.d_*pr.d_*pr.w_;
	if (max_ < pr.d_) { max_ = pr.d_; nmax_ = i; }
}

void mcHistogram::operator+=(const mcHistogram& h)
{
	if (comulative_)
		throw std::exception("попытка добавления данных к интегральной гистограмме");
	if (h.npnt() != npnt() || h.argmin() != argmin() || h.step() != step())
		throw std::exception("несогласованные размеры гистограмм");
	for (unsigned i = 0; i < npnt(); i++) {
		d_[i] += h.d()[i];
		total_ += h.d()[i];
	}
	sum_ += h.sum();
	sum2_ += h.sum2();
	max_ = (max_ < h.argmax()) ? h.argmax() : max_;
	nmax_ = (nmax_ < h.nmax()) ? h.nmax() : nmax_;
}

ostream& operator << (ostream& os, const mcHistogram& h)
{
	os << "NP =" << '\t' << h.nmax() << endl;
	for (unsigned i = 0; i < h.npnt(); i++)
		os << h.arg(i) << '\t' << h.d()[i] << endl;
	return os;
}

ostream& operator << (ostream& os, const vector<mcHistogram>& ha)
{
	unsigned int i, j, nh = (int)ha.size();
	if (!nh) return os;  // skip empty
	unsigned int np = ha[0].npnt();
	if (!np) return os;  // skip empty
	os << "nhistograms = " << nh << '\t' << "np =" << '\t' << np << endl;
	for (j = 0; j < np; j++) {
		os << ha[0].arg(j);
		for (i = 0; i < nh; i++) {
			if (ha[i].npnt() != np) throw std::exception("несогласованные размеры гистограмм");
			os << '\t' << ha[i].d()[j];
		}
		os << endl;
	}
	return os;
}
