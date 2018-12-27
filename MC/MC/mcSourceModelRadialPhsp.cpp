#include "mcSourceModelRadialPhsp.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "mcPhysics.h"
#include <math.h>

mcSourceModelRadialPhsp::mcSourceModelRadialPhsp(const char* name, int nThreads, double z0)
:mcSource(name, nThreads), isxray_(0), z0_(z0), ne_(0), nr_(0), nthet_(0), naxial_(0), 
focus_(0), emax_(0), rmax_(0), thetmax_abs_(0), thetmax_(0), 
de_(0), dr_(0), dthet_abs_(0), dthet_(0), daxial_(0), 
we117_(0), we133_(0), wexscale_(1.0)
{
}

void mcSourceModelRadialPhsp::sample(mcParticle& p, mcThread* thread)
{
	sample(p, thread->rng());

	p.plast = p.p;
	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke;
}

void mcSourceModelRadialPhsp::sample(mcParticle& p, mcRng& rng) const
{
	p.t = MCP_PHOTON;
	p.q = 0;

	// Энергия
	double rnd = rng.rnd();
	double rnd0 = rng.rnd();
	double rnd1 = rng.rnd();
	double rnd2 = rng.rnd();
	double f;	// коэффициент интерполяции
	double r, t, a;
	bool isGammaPrimary = true;;

	if(rnd < we117_ && isxray_ == 0)
	{
		p.ke = 1.17;
		sampleForEnergy(ne_, r, t, a, rnd0, rnd1, rnd2);
	}
	else if (rnd < we133_ && isxray_ == 0)
	{
		p.ke = 1.33;
		sampleForEnergy(ne_ + 1, r, t, a, rnd0, rnd1, rnd2);
	}
	else
	{
		isGammaPrimary = false;;
		p.ke = xehistogram_.sample((1.0 - rnd) * wexscale_);
		int eidx = int(p.ke / de_);
		if(eidx >= ne_ - 1) eidx = ne_ - 2;
		int eidx2 = eidx + 1;
		f = p.ke / de_ - double(eidx);

		// GG 20160212 После предположения о том, что проблемы могут создаваться
		// при интерполяции по родительским измерениям из-за больших шумов 
		// решено исключить интерполяцию и брать параметр для ближайшего родителя.

		//double r1, t1, a1;
		//double r2, t2, a2;
		//sampleForEnergy(eidx, r1, t1, a1, rnd0, rnd1, rnd2);
		//sampleForEnergy(eidx2, r2, t2, a2, rnd0, rnd1, rnd2);
		//r = (1 - f) * r1 + f * r2;
		//t = (1 - f) * t1 + f * t2;
		//a = (1 - f) * a1 + f * a2;
		sampleForEnergy(f < 0.5 ? eidx : eidx2, r, t, a, rnd0, rnd1, rnd2);
	}

	p.p.set(-r, 0, 0);

	geomVector3D vp(-r, 0, focus_);
	vp.normalize();
	double cost = isGammaPrimary ? 1 - t : cos(t);
	double sint = sqrt(1 - cost * cost);
	double cosa = cos(a);
	double sina = sin(a);
	if(rng.rnd() < 0.5) sina = -sina;

	p.u.set(cost*vp.x() + sint*cosa * vp.z(), sint*sina, cost*vp.z() - sint*cosa * vp.x());

	// Поворачиваем на случайный угол вокруг оси системы.
	double cosPhi, sinPhi;
	mcPhysics::GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);

	p.p.p_[1] = -p.p.p_[0] * sinPhi;
	p.p.p_[0] *= cosPhi;
	p.p.p_[2] += z0_;

	double x = p.u.p_[0], y = p.u.p_[1];
	p.u.p_[0] = x * cosPhi + y * sinPhi;
	p.u.p_[1] = -x * sinPhi + y * cosPhi;
}

void mcSourceModelRadialPhsp::sampleForEnergy(int eidx, double& r, double& t, double& a, double rnd0, double rnd1, double rnd2) const
{
	double f = rhistograms_[eidx].sample(rnd0);
	int ridx = int(f / dr_), ridx2 = ridx + 1;

	// Экстраполяция может приводить к фатальным проблеммам
	if (ridx2 >= nr_) ridx2 = nr_ - 1;
	double fr = f / dr_ - double(ridx);
	r = rmax_ * acos(1 - f);

	// GG 20160212 Отказ от интерполяции.

	//double t1 = thistograms_[eidx][ridx].sample(rnd1);
	//double t2 = thistograms_[eidx][ridx2].sample(rnd1);
	//t = (1 - fr) * t1 + fr * t2;
	t = thistograms_[eidx][fr < 0.5 ? ridx : ridx2].sample(rnd1);

	double dthet = eidx < ne_ ? dthet_abs_ : dthet_;
	int tidx = int(t / dthet), tidx2 = tidx + 1;
	if (tidx2 >= nthet_) tidx2 = nthet_ - 1;
	double ft = t / dthet - double(tidx);

	if (eidx >= ne_)
	{
		t = thetmax_ * asin(t) * 0.25;
	}
	//else
	//{
	//	t = thetmax_abs_ * t;
	//}

#ifdef SRC_USE_ASIMUT_MODEL
	//double a11 = sampleAzimut(apars_[eidx][ridx][tidx], rnd2);
	//double a12 = sampleAzimut(apars_[eidx][ridx][tidx2], rnd2);
	//double a21 = sampleAzimut(apars_[eidx][ridx2][tidx], rnd2);
	//double a22 = sampleAzimut(apars_[eidx][ridx2][tidx2], rnd2);

	a = sampleAzimut(apars_[eidx][fr < 0.5 ? ridx : ridx2][ft < 0.5 ? tidx : tidx2], rnd2);
#else
	//double a11 = ahistograms_[eidx][ridx][tidx].sample(rnd2);
	//double a12 = ahistograms_[eidx][ridx][tidx2].sample(rnd2);
	//double a21 = ahistograms_[eidx][ridx2][tidx].sample(rnd2);
	//double a22 = ahistograms_[eidx][ridx2][tidx2].sample(rnd2);

	a = ahistograms_[eidx][fr < 0.5 ? ridx : ridx2][ft < 0.5 ? tidx : tidx2].sample(rnd2);
#endif

	//a = (1 - fr)*((1 - ft)*a11 + ft*a12) + fr*((1 - ft)*a21 + ft*a22);
}

double mcSourceModelRadialPhsp::sampleAzimut(const vector<double>& a, double x)
{
	return x * (a[0] * (x*x - 1) + a[1] * (x - 1) + PI);
	//return a[0] == 0 ? PI * x : (sqrt(1.0 + a[0] * x) - 1.0) * a[1];
}

void* mcSourceModelRadialPhsp::saveToMemory(int& size)
{
	// Размер необходимой памяти
	int memcount = 0;
	memcount += sizeof(int);		// тип излучения
	memcount += 4 * sizeof(int);	// размерность 4D пространства
	memcount += 5 * sizeof(double);	// фокус и максимальные значения координат (кроме фиксированного азимута)
	memcount += 3 * sizeof(double);	// веса монолиний и рассеянного спектра излучения

	int ne = isxray_ == 0 ? ne_ + 2 : ne_;
	memcount += (ne_ + 1) * sizeof(double);							// гистограмма спектра рассеянного излучения
	memcount += ne * (nr_ + 1) * sizeof(double);					// гистограммы распределений по радиусу
	memcount += ne * nr_ * (nthet_ + 1) * sizeof(double);			// гистограммы распределений по углу
#ifdef SRC_USE_ASIMUT_MODEL
	memcount += ne * nr_ * nthet_ * 2 * sizeof(double);	// гистограммы распределений по азимуту
#else
	memcount += ne * nr_ * nthet_ * (naxial_ + 1) * sizeof(double);	// гистограммы распределений по азимуту
#endif

	void* buffer = malloc(memcount);
	if(buffer == nullptr)
		throw std::exception("mcSourceModelRadialPhsp::saveToMemory can not allocate memory");

	char* pbuffer = (char*)buffer;
	*(int*) pbuffer = isxray_; pbuffer += sizeof(int);
	*(int*) pbuffer = ne_; pbuffer += sizeof(int);
	*(int*) pbuffer = nr_; pbuffer += sizeof(int);
	*(int*)pbuffer = nthet_; pbuffer += sizeof(int);
	*(int*)pbuffer = naxial_; pbuffer += sizeof(int);

	*(double*)pbuffer = focus_; pbuffer += sizeof(double);
	*(double*)pbuffer = emax_; pbuffer += sizeof(double);
	*(double*)pbuffer = rmax_; pbuffer += sizeof(double);
	*(double*)pbuffer = thetmax_abs_; pbuffer += sizeof(double);
	*(double*)pbuffer = thetmax_; pbuffer += sizeof(double);

	*(double*)pbuffer = we117_; pbuffer += sizeof(double);
	*(double*)pbuffer = we133_; pbuffer += sizeof(double);
	*(double*)pbuffer = wexscale_; pbuffer += sizeof(double);

	int len = xehistogram_.npnt() * sizeof(double);
	memcpy(pbuffer, xehistogram_.data(), len);
	pbuffer += len;

	int ie;
	for(ie = 0; ie<ne; ie++)
	{
		len = rhistograms_[ie].npnt() * sizeof(double);
		memcpy(pbuffer, rhistograms_[ie].data(), len);
		pbuffer += len;
	}

	for(ie = 0; ie<ne; ie++)
	{
		for(int ir = 0; ir<nr_; ir++)
		{
			len = thistograms_[ie][ir].npnt() * sizeof(double);
			memcpy(pbuffer, thistograms_[ie][ir].data(), len);
			pbuffer += len;
		}
	}

	for(ie = 0; ie<ne; ie++)
	{
		for(int ir = 0; ir<nr_; ir++)
		{
			for(int it = 0; it<nthet_; it++)
			{
#ifdef SRC_USE_ASIMUT_MODEL
				*(double*)pbuffer = apars_[ie][ir][it][0]; pbuffer += sizeof(double);
				*(double*)pbuffer = apars_[ie][ir][it][1]; pbuffer += sizeof(double);
#else
				len = ahistograms_[ie][ir][it].npnt() * sizeof(double);
				memcpy(pbuffer, ahistograms_[ie][ir][it].data(), len);
				pbuffer += len;
#endif
			}
		}
	}

	size = memcount;
	return buffer;
}

void mcSourceModelRadialPhsp::readFromMemory(void* buffer)
{
	char* p = (char*)buffer;
	isxray_ = *(int*) (p); p += sizeof(int);
	ne_ = *(int*) (p); p += sizeof(int);
	nr_ = *(int*) (p); p += sizeof(int);
	nthet_ = *(int*)(p); p += sizeof(int);
	naxial_ = *(int*)(p); p += sizeof(int);

	focus_ = *(double*)(p); p += sizeof(double);
	emax_ = *(double*)(p); p += sizeof(double);
	rmax_ = *(double*)(p); p += sizeof(double);
	thetmax_abs_ = *(double*) (p); p += sizeof(double);
	thetmax_ = *(double*)(p); p += sizeof(double);

	we117_ = *(double*)(p); p += sizeof(double);
	we133_ = *(double*)(p); p += sizeof(double);
	wexscale_ = *(double*)(p); p += sizeof(double);

	de_ = emax_ / ne_;
	dr_ = 2.0 / nr_;
	dthet_abs_ = thetmax_abs_ / nthet_;
	dthet_ = 1.0 / nthet_;
	daxial_ = PI / naxial_;

	createHistograms(ne_, nr_, nthet_, naxial_);

	xehistogram_.restore(ne_ + 1, (double*)(p)); p += (ne_ + 1) * sizeof(double);

	int ie, ne = isxray_ == 0 ? ne_ + 2 : ne_;
	for(ie = 0; ie<ne; ie++)
	{ rhistograms_[ie].restore(nr_ + 1, (double*)(p)); p += (nr_ + 1) * sizeof(double); }

	for(ie = 0; ie<ne; ie++)
	{
		for(int ir = 0; ir<nr_; ir++)
		{ 
			thistograms_[ie][ir].restore(nthet_ + 1, (double*)(p)); p += (nthet_ + 1) * sizeof(double); 
		}
	}

	for(ie = 0; ie<ne; ie++)
	{
		for(int ir = 0; ir<nr_; ir++)
		{
			for(int it = 0; it<nthet_; it++)
			{ 
#ifdef SRC_USE_ASIMUT_MODEL
				apars_[ie][ir][it][0] = *(double*)p; p += sizeof(double); 
				apars_[ie][ir][it][1] = *(double*)p; p += sizeof(double); 
#else
				ahistograms_[ie][ir][it].restore(naxial_ + 1, (double*)(p)); p += (naxial_ + 1) * sizeof(double);
#endif
			}
		}
	}
}

void mcSourceModelRadialPhsp::createHistograms(int ne, int nr, int nt, int na)
{
	ne_ = ne;
	nr_ = nr;
	nthet_ = nt;
	naxial_ = na;

	if(isxray_ == 0)
		ne += 2;	// Поддержка кобальта с двумя монолиниями

	rhistograms_.resize(ne);
	thistograms_.resize(ne);

#ifdef SRC_USE_ASIMUT_MODEL
	apars_.resize(ne);
#else
	ahistograms_.resize(ne);
#endif
	for(int ie = 0; ie < ne; ie++)
	{
		thistograms_[ie].resize(nr_);
#ifdef SRC_USE_ASIMUT_MODEL
		apars_[ie].resize(nr_);
#else
		ahistograms_[ie].resize(nr_);
#endif
		for(int ir = 0; ir < nr_; ir++)
		{
#ifdef SRC_USE_ASIMUT_MODEL
			apars_[ie][ir].resize(nthet_);
			for(int it = 0; it < nthet_; it++)
				apars_[ie][ir][it].resize(2);
#else
			ahistograms_[ie][ir].resize(nthet_);
#endif
		}
	}
}

void mcSourceModelRadialPhsp::createFromDistribution(int isxray, int ne, int nr, int nt, int na,
													 double emax, double rmax, double tmax, double focus,
													 const double* data)
{
	isxray_ = isxray;
	createHistograms(ne, nr, nt, na);

	emax_ = emax;
	de_ = emax_ / ne_;
	rmax_ = rmax / PI;
	thetmax_abs_ = tmax * PI / 180.0;
	thetmax_ = 2.0 * (1 - cos(thetmax_abs_)) / PI;
	focus_ = focus;

	de_ = emax_ / ne_;
	dr_ = 2.0 / nr_;
	dthet_abs_ = thetmax_abs_ / nthet_;
	dthet_ = 1.0 / nthet_;
	daxial_ = PI / naxial_;

	if (isxray_ == 0)
		ne += 2;	// Поддержка кобальта с двумя монолиниями

	int nnt = nt * na;
	int nnr = nr * nnt;

	std::vector<double> edistribution(ne, 0);
	double ecount = 0;

	for(int ie = 0; ie < ne; ie++)
	{
		double rcount = 0;
		std::vector<double> rdistribution(nr, 0);

		for(int ir = 0; ir < nr; ir++)
		{
			std::vector<double> tdistribution(nthet_, 0);
			double tcount = 0;

			for(int it = 0; it < nt; it++)
			{
				double acount = 0;
				for(int ia = 0; ia < na; ia++)
					acount += data[ie * nnr +  ir * nnt + it * na + ia];

#ifdef SRC_USE_ASIMUT_MODEL
				double a, b;
				fitAsimutDistribution(na, data + ie * nnr +  ir * nnt + it * na, a, b);
				apars_[ie][ir][it][0] = a;
				apars_[ie][ir][it][1] = b;
#else
				ahistograms_[ie][ir][it].setFromDistribution(na, na, 0, daxial_, data + ie * nnr +  ir * nnt + it * na);
#endif
				tdistribution[it] = acount;
				tcount += acount;
			}

			rdistribution[ir] = tcount;
			rcount += tcount;
			thistograms_[ie][ir].setFromDistribution(nt, 0, ie < ne_ ? dthet_abs_ : dthet_, tdistribution);
		}

		rhistograms_[ie].setFromDistribution(nr, 0, dr_, rdistribution);
		if(ie == ne_)
			we117_ = rcount;
		else if(ie == ne_ + 1)
			we133_ = rcount;
		else
			edistribution[ie] = rcount;
		ecount += rcount;
	}

	xehistogram_.setFromDistribution(ne_, 0, de_, edistribution);
	if (isxray_ == 0)
	{
		we117_ /= ecount;
		we133_ /= ecount;
		we133_ += we117_;
		wexscale_ = 1.0 / (1.0 - we133_);
	}
}

void mcSourceModelRadialPhsp::fitAsimutDistribution(int np, const double* data, double& a, double& b)
{
	// GG 20160216 Это старая первая попытка апроксимировать распределение по азимуту аналитически
	//--------------------------------------------------------------------------------------------
	// Среднестатистические значения для рассеянных фотнов
	a = 1.6539;
	b = -2.9756;

	int i;
	vector<double> c(np);
	double sum = 0;
	c[0] = data[0];
	for(i=1; i < np; i++)
		c[i] = c[i-1] + data[i];

	if(c[np-1] > 0)
	{
		double f = 1.0 / c[np-1];
		for(i=0; i < np; i++)
			c[i] *= f;

		// Подгоняем методом наименьших квадратов функцию
		// f(x) = x * [a * (x^2-1) + b * (x-1) + PI]
		// Если плотность распределения апроксимировать параборлой,
		// то в комулятивном распределении, представлямом интегралом от плотности,
		// появляется первый множитель X.
		// Если еще потребовать прохождение через точку PI (диапазон косинусов) / 1,
		// то один коэффициент исчезает и появляется последняя двойка.

		// Коэффициенты системы двух линейных уравнений
		double A1 = 0, AB = 0, B2 = 0, F1 = 0, F2 = 0;
		double danorm = PI / np;

		for(i=0; i<np; i++)
		{
			double x = c[i];
			double x1 = x - 1;
			double x21 = x * x - 1;
			f = danorm * (i + 1);

			A1 += x * x * x21 * x21;
			AB += x * x * x1 * x21;
			B2 += x * x * x1 * x1;
			F1 += x * x21 * (f - PI * x);
			F2 += x * x1 * (f - PI * x);
		}

		double C = A1 * B2 - AB * AB;
		if(C != 0)
		{
			a = (F1 * B2 - F2 * AB) / C;
			b = (F2 * A1 - F1 * AB) / C;
		}
	}
	/*
	*/

	// GG 20160216 Это новая попытка на основе прямой апроксимации 
	// обратнойкомулятивной функции распределения
	//-------------------------------------------------------------

	/*
	// Значения по умолчанию соответствуют равномерному распределению во всем диапазоне.
	// Это определяется по нулевому значению параметра a в самплинге.
	a = b = 0;

	// Конвертируем распределение в комулятивное и нормируем его на 1 по обеим координатам. 
	int i;
	vector<double> c(np + 1);
	c[0] = 0;
	for (i = 0; i < np; i++)
		c[i+1] = c[i] + data[i];

	if (c[np] > 0)
	{
		// На входе мы получаем только количество частиц в бине.
		// Массив разбит на np точек с постоянным шагом.
		// Особенность в том, что для простоты рассмотрение мы предполагаем весь диапазон равен не Пи, а 1.
		double danorm = 1.0 / np;

		double f = 1.0 / c[np];
		for (i = 0; i <= np; i++)
			c[i] *= f;
		
		double A = 0, B = 0;

		// Точки 0 и 1 опускаем, так как их роложение предопределено самой формулой апроксимации
		// (две неизвестные мы превратили в одну, заставив апроксимацию проходить через них)
		for (i = 1; i<np; i++)
		{
			double x = i * danorm;
			double y = c[i];
			A += (y - x) * x * (x - 1);
			double bi = x * (x - 1);
			B += bi * bi;
		}

		if (B != 0)
		{
			f = A / B;
			// При малой статистике возможны отрицательные значения f.
			// Закрываем на это глаза считая распределение равномерным.
			if (f > 0.0 && f < 1.0)
			{
				a = 4 * f / ((1 - f) * (1 - f));
				b = PI * (1 - f) / (2 * f);
			}
		}
	}
	*/
}

ostream& operator << (ostream& os, const mcSourceModelRadialPhsp& s)
{
	os << (const mcSource&)s;
	os << "TYPE = \t" << s.isxray_ << endl;
	os << "NAME = \t" << s.getName() << endl;
	os << "POSITION = \t" << s.z0_ << endl;
	os << endl;
	os << "z0_ = \t" << s.z0_ << endl;
	os << "ne_ = \t" << s.ne_ << endl;
	os << "nr_ = \t" << s.nr_ << endl;
	os << "nthet_ = \t" << s.nthet_ << endl;
	os << "naxial_ = \t" << s.naxial_ << endl;
	os << endl;
	os << "focus_ = \t" << s.focus_ << endl;
	os << "emax_ = \t" << s.emax_ << endl;
	os << "rmax_ = \t" << s.rmax_ * PI << endl;
	os << "thetmax_abs_ = \t" << s.thetmax_abs_ * 180 / PI << endl;
	os << endl;
	os << "de_ = \t" << s.de_ << endl;
	os << "dr_ = \t" << s.dr_ << endl;
	os << "dthet_abs_ = \t" << s.dthet_abs_ << endl;
	os << "dthet_ = \t" << s.dthet_ << endl;
	os << "daxial_ = \t" << s.daxial_ << endl;
	os << endl;
	os << "we117_ = \t" << s.we117_ << endl;
	os << "we133_ = \t" << s.we133_ << endl;
	os << "wexscale_ = \t" << s.wexscale_ << endl;

	os << endl;
	os << "Распределение по энергиям:" << endl;
	os << "--------------------------" << endl;
	os << s.xehistogram_;

	os << endl;
	os << "Распределения по радиусам:" << endl;
	os << "--------------------------" << endl;
	os << s.rhistograms_;

	os << endl;
	os << "Распределения по углам:" << endl;
	os << "--------------------------" << endl;
	int i;
	for (i = 0; i < (s.isxray_ == 0 ? s.ne_ + 2 : s.ne_); i++)
	{
		os << endl << "EBIN = \t" << i << endl << endl;
		os << s.thistograms_[i];
	}
	
#ifdef SRC_USE_ASIMUT_MODEL
	int j, k;
	os << endl;
	os << "Распределения по азимутам:" << endl;
	for(i = 0; i < (s.isxray_ == 0 ? s.ne_ + 2 : s.ne_); i++)
	{
		os << endl << "EBIN = \t" << i << endl << endl;
		os << "--------------------------" << endl;
		os << "Коэффициенты a:" << endl;
		os << "--------------------------" << endl;
		for(k = 0; k < s.nthet_; k++)
			os << "\t" << (k+0.5) * s.dthet_;
		os << endl;
		for(j = 0; j < s.nr_; j++)
		{
			os << (j + 0.5) * s.dr_;
			for(k = 0; k < s.nthet_; k++)
				os << "\t" << s.apars_[i][j][k][0];
			os << endl;
		}

		os << endl;
		os << "Коэффициенты b:" << endl;
		os << "--------------------------" << endl;
		for(k = 0; k < s.nthet_; k++)
			os << "\t" << (k+0.5) * s.dthet_;
		os << endl;
		for(j = 0; j < s.nr_; j++)
		{
			os << (j + 0.5) * s.dr_;
			for(k = 0; k < s.nthet_; k++)
				os << "\t" << s.apars_[i][j][k][1];
			os << endl;
		}
	}
#endif

	return os;
}

void mcSourceModelRadialPhsp::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	//Источник задается в абсолютных координатах.
	//Поэтому ему не нужна матрица преобразования координат.

	int ir, it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;

	os << "# Source: " << name_ << endl;

	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	// Концентрические круги
	for(ir=1; ir<=nr_; ir++) {
		double r = rmax_ * acos(1 - dr_ * ir);
		for(it=0; it<360; it+=da) {
			geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), z0_);
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r*sin(mPi*(it+da)), r*cos(mPi*(it+da)), z0_);
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			count++;
		}
	}

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	for(it=0; it<count; it++)
		os << "      " << 2*it << ' ' << 2*it+1 << " -1" << endl;
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;
}
