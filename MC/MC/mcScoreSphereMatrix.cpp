#include "mcScoreSphereMatrix.h"
#include "mcTransport.h"

mcScoreSphereMatrix::mcScoreSphereMatrix(const char* module_name, int nThreads, int np, int nm, double rmin, double rmax)
	:mcScore(module_name, nThreads),
	np_(np), nm_(nm), rmin_(rmin), rmax_(rmax)
{
	int len = nThreads * np_ * nm_;
	MAll_ = new double[len];
	memset(MAll_, 0, len * sizeof(double));
	M_ = new double*[nThreads];
	M_[0] = MAll_;
	for (int i = 1; i < nThreads_; i++)
		M_[i] = M_[i - 1] + np_ * nm_;

	r2min_ = rmin_ * rmin_;
	r2max_ = rmax_ * rmax_;
	pstep_ = 2. / np_;	// мы разбиваем косину (p.z()) на равные интервалы
	mstep_ = 2. * PI / nm_;
}

mcScoreSphereMatrix::~mcScoreSphereMatrix()
{
	delete[] MAll_;
	delete[] M_;
}

void mcScoreSphereMatrix::ScoreFluence(const mcParticle& particle) { }

void mcScoreSphereMatrix::ScorePoint(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p)
{
	double r2 = p.sqLength();
	if (r2 < r2min_ || r2 > r2max_) return;

	geomVector3D pp(p);
	pp.normalize();
	int ip = int((1.0 - pp.z()) / pstep_);
	if (ip < 0 || ip >= np_) return;

	double xy = p.lengthXY();
	if (xy == 0) return;	// сверхмаловероятно, просто защищаем код
	int im = int(acos(p.x()/xy) / mstep_);
	if (p.y() < 0) im = nm_ - im - 1;
	if (im < 0 || im >= nm_) return;

	scoreEnergyInVoxel(iThread, ip, im, edep);
}

void mcScoreSphereMatrix::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)								   
{
	// Если обе точки внутри внутренней сферы, то нет контакта с детектором.
	double r0 = p0.length();
	double r1 = p1.length();
	if (r0 <= rmin_ && r1 <= rmin_) return;

	double length = (p1 - p0).length();
	double eddens = edep / length;

	// Возможна ситуация, когда два отрезка данного отрезка лежат внутри сферического кольца.
	// Такую ситуацию можно идентифицировать по расстоянию до ближайшей точки отрезка.
	// В случае двоения формируем два отрезка и вызываем эту функцию для каждого рекурентно.

	// Сортируем точки по возрастанию радиуса
	geomVector3D ps0(p0), ps1(p1);	// отсортированные точки
	if (r0 > r1)
	{
		ps0 = p1;
		ps1 = p0;
		double f = r1; r1 = r0; r0 = f;
	}

	// Q: Как определить минимальное расстояние от центра координат до отрезка в пространстве?
	// Строим вектор нормали из центра к прямой линии
	geomVector3D p3 = ps0 ^ (ps1 - ps0);
	geomVector3D n = (ps1 - ps0) ^ p3;
	n.normalize();
	double h = ps0 * n;
	if (h >= rmax_) return;

	// Единичный вектор направления отрезка
	geomVector3D p01 = ps1 - ps0;
	p01.normalize();

	// Если высота между радиусами до точек отрезка, то деление возможно
	//if (h > r0 && h < r1)
	if ((p01 * ps0) * (p01 * ps1) <0)
	{
		// Деление внешним радиусом. Остается один отрезок, но короче.
		// Обрезаем и продолжаем распределение.
		if (r0 > rmax_ && r1 > rmax_)
		{
			double hp = sqrt(r2max_ - h*h);
			double x = sqrt(r0*r0 - h*h) - hp;
			length -= x;
			ps0 += p01 * x;

			x = hp - sqrt(r1*r1 - h*h);
			length += x;
			ps1 += p01 * x;

			r0 = r1 = rmax_;
		}
		
		// Деление внутренним радиусом
		if (h < rmin_ && r0 > rmin_ && r1 > rmin_)
		{
			double hp = sqrt(r2min_ - h*h);
			double x = sqrt(r0*r0 - h*h) - hp;
			geomVector3D pend = ps0 + (p01 * x);
			ScoreLine(eddens * x, iThread, region, pt, ps0, pend);

			x = sqrt(r1*r1 - h*h) - hp;
			geomVector3D pstart = ps1 - (p01 * x);
			ScoreLine(eddens * x, iThread, region, pt, pstart, ps1);

			return;
		}
	}
	if (r0 > rmax_) return;
		
	// Проверяем и при необходимости обрезаем в случае однократного пересечения	какой либо сферы.
	// Мы точно знаем, что радиусы отсортированы по возрастанию и возможно не более одного пересечения.
	if (r1 > rmax_)
	{
		double hp = sqrt(r2max_ - h*h);
		double x = hp - sqrt(r1*r1 - h*h);
		length += x;
		ps1 += p01 * x;
		r1 = rmax_;
	}
	if (r0 < rmin_)
	{
		double hp = sqrt(r2min_ - h*h);
		double x = sqrt(r0*r0 - h*h) - hp;
		length += x;
		ps0 += p01 * x;
		r0 = rmin_;
	}

	// Теперь мы гарантированно имеем отрезок между ограничивающими сферами.
	// Приступаем к распределению между ячейками двумерной матрицы.

	// Индексы начала
	geomVector3D pps0(ps0);
	pps0.normalize();
	int ip0 = int((1.0 - pps0.z()) / pstep_);
	if (ip0 < 0 || ip0 >= np_) return;			// Hack!! Возможная защита кода

	double xy = ps0.lengthXY();
	if (xy == 0) return;	// сверхмаловероятно, просто защищаем код
	int im0 = int(acos(ps0.x() / xy) / mstep_);
	if (ps0.y() < 0) im0 = nm_ - im0 - 1;
	if (im0 < 0 || im0 >= nm_) return;

	bool isMPositive = (ps0 ^ ps1).z() > 0;	// направление отрезка по меридианам
	bool isPPositive = false;
	bool willParallelCross = true;
	bool willMeridianCross = true;
	bool needParallelCross = true;
	bool needMeridianCross = true;
	double dp = 0, dm = 0;

	// Идем по отрезку отслеживая границы вокселей
	while (true)
	{
		if (willParallelCross && needParallelCross)
		{
			double dp1, dp2;
			// Индикатор в каком полушарии находимся
			double cosz = 1. - (ip0 + 1) * pstep_;
			if (cosz >= 0)
			{
				dp1 = distanceToParallel(ps0, p01, ip0, true);
				dp2 = distanceToParallel(ps0, p01, ip0 + 1, false);
			}
			else
			{
				dp1 = distanceToParallel(ps0, p01, ip0, false);
				dp2 = distanceToParallel(ps0, p01, ip0 + 1, true);
			}
			isPPositive = dp2 < dp1;
			dp = fmin(dp1, dp2);
			if (dp >= length)
				willParallelCross = false;
		}
		else
			needParallelCross = false;

		if (willMeridianCross && needMeridianCross)
		{
			int idx = isMPositive ? (im0 + 1) % nm_ : im0;
			dm = distanceToMeridian(ps0, p01, idx, isMPositive);
			if (dm >= length)
				willMeridianCross = false;
		}
		else
			needMeridianCross = false;

		// Не дотягиваем ни до одной границы
		if (!willParallelCross && !willMeridianCross)
		{
			scoreEnergyInVoxel(iThread, ip0, im0, eddens * length);
			break;
		}

		if (dm > dp)
		{
			scoreEnergyInVoxel(iThread, ip0, im0, eddens * dp);
			ps0 += p01 * dp;
			ip0 = isPPositive ? (ip0 + 1) : (ip0 - 1);
			length -= dp;
			dm -= dp;
			needParallelCross = true;
		}
		else
		{
			scoreEnergyInVoxel(iThread, ip0, im0, eddens * dm);
			ps0 += p01 * dm;
			im0 = (isMPositive ? (im0 + 1) : (nm_ + im0 - 1)) % nm_;
			length -= dm;
			dp -= dm;
			needMeridianCross = true;
		}
	}
}

double mcScoreSphereMatrix::distanceToParallel(const geomVector3D& p, const geomVector3D& v, int ip, bool isPOutside)
{
	// Полюса, с которыми нельзя столкнуться
	if (ip <= 0 || ip >= np_)
		return DBL_MAX;

	// Определение наклона конуса и полушария
	double cosz = 1. - ip * pstep_;
	double cos2z = cosz * cosz;

	// Решение квадратного уравнения (ad^2 + 2bd + c = 0).
	// Коэффициенты уравнения не зависят от того, откуда летит частица,
	// но правильное решение зависит от множества условий.
	double vx = v.x(), vy = v.y(), vz = v.z();
	double px = p.x(), py = p.y(), pz = p.z();
	double a = vz * vz - cos2z;
	double b = vz * pz - cos2z * (vx * px + vy * py + vz * pz);
	double c = pz * pz - cos2z * (px * px + py * py + pz * pz);
	
	double det = b * b - a * c;
	if (det < 0)	// летим мимо
		return DBL_MAX;

	if(isPOutside)
	{
		// Если пересечений нет вообще, то ситуация обработана выше через детерминант.
		// Возможно одно пересечение
		if (fabs(a) <= FLT_EPSILON)
		{
			if (b == 0) return DBL_MAX;
			double d = -0.5 * c / b;
			if(d <= 0) return DBL_MAX;
			else return d;
		}
		else
		{
			double d1 = (-b - sqrt(det)) / a;
			double d2 = (-b + sqrt(det)) / a;
			if(d1 <= 0 && d2 <= 0) return DBL_MAX;
			double d = (d1 > 0 && d2 > 0) ? fmin(d1, d2) : d1 > 0 ? d1 : d2;
			// Последний фильтр - проверка  на пересечение правильного полушария
			double f = (p + v * d).z();
			if (f * cosz > 0)
				return d;
			else
				return DBL_MAX;
		}
	}
	else
	{
		double d = a == 0 ? (b == 0 ? -1 : -0.5 * c / b) : (-b - sqrt(det)) / a;
		if (d <= 0)
			return DBL_MAX;	// решение есть, но оно не в той стороне куда летим
		else return d;
	}
}

double mcScoreSphereMatrix::distanceToMeridian(const geomVector3D& p0, const geomVector3D& v, int im, bool isMPositive)
{
	// Нормаль к плоскости меридиана
	geomVector3D n(-sin(im * mstep_), cos(im * mstep_), 0);
	double f = n * v;
	if ((isMPositive && f <= 0) || (!isMPositive && f >= 0))
		return DBL_MAX;
	else
		return -(p0 * n) / f;
}

void mcScoreSphereMatrix::scoreEnergyInVoxel(int iThread, int ip, int im, double edep)
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreSphereMatrix::scoreEnergyInVoxel: thread exceed container size");

	if (ip >= 0 && ip < np_ && im >= 0 && im < nm_) 
	{
		M_[iThread][ip*nm_ + im] += edep;
		etotal_[iThread] += edep;
	}
}

void mcScoreSphereMatrix::CE2D()
{
	if (dconverted_)
		throw std::exception("mcScoreSphereMatrix::CE2D: Already converted to dose");
	dconverted_ = true;

	double dr3 = (rmax_*r2max_ - rmin_*r2min_) / 3.0;
	double vol = mstep_ * pstep_ * dr3;

	for (int ip = 0; ip < np_; ip++)
	{
		for (int im = 0; im < nm_; im++)
		{
			for (int i = 0; i < nThreads_; i++)
				M_[i][ip*nm_ + im] /= vol;
		}
	}
}

double mcScoreSphereMatrix::Dose(int ip, int im) const
{
	double f = 0;
	if (ip >= 0 && ip < np_ && im >= 0 && im < nm_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += M_[i][ip*nm_ + im];
	}
	return f;
}

void mcScoreSphereMatrix::dumpVRML(ostream& os) const
{
	os << "# mcScoreSphereMatrix Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int count = 0;
	double mPi = PI / 180;

	os << "# mcScoreConicalRZ Score: " << name_ << endl;
	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	for (int ip = 0; ip < np_; ip++)
	{
		double m = 1 - (ip * pstep_);
		double z1 = rmin_ * m, z2 = rmax_ * m;
		double z3 = rmin_ * (m - pstep_), z4 = rmax_ * (m - pstep_);
		double r1 = sqrt(rmin_ * rmin_ - z1 * z1);
		double r2 = sqrt(rmax_ * rmax_ - z2 * z2);
		double r3 = sqrt(rmin_ * rmin_ - z3 * z3);
		double r4 = sqrt(rmax_ * rmax_ - z4 * z4);

		for (int im = 0; im < nm_; im++)
		{
			double a1 = im * mstep_;
			double a2 = (im + 1) * mstep_;

			geomVector3D p = geomVector3D(r1 * sin(a1), r1 * cos(a1), z1) * mttow;	// 0
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r1 * sin(a2), r1 * cos(a2), z1) * mttow;				// 1
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r3 * sin(a1), r3 * cos(a1), z3) * mttow;				// 2
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r2 * sin(a1), r2 * cos(a1), z2) * mttow;				// 3
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r2 * sin(a2), r2 * cos(a2), z2) * mttow;				// 4
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r4 * sin(a1), r4 * cos(a1), z4) * mttow;				// 5
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			
			count += 6;
		}
	}

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	for (int it = 0; it < count; it += 6)
	{
		os << "      " << it << ' ' << it + 1 << " -1" << endl;
		os << "      " << it << ' ' << it + 2 << " -1" << endl;
		os << "      " << it << ' ' << it + 3 << " -1" << endl;
		os << "      " << it + 3 << ' ' << it + 4 << " -1" << endl;
		os << "      " << it + 3 << ' ' << it + 5 << " -1" << endl;
	}
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;

}

void mcScoreSphereMatrix::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	os << "Spherical dose matrix:" << endl;
	os << endl;

	int ip, im;
	for (im = 0; im < nm_; im++)
		os << "\t" << ((im + 0.5) * mstep_ * 180. / PI);
	os << endl;
	
	for (ip = 0; ip < np_; ip++)
	{
		os << (acos(1.0 - (ip + 0.5) * pstep_) * 180 / PI);
		for (im = 0; im < nm_; im++)
			os << "\t" << Dose(ip, im);
		os << endl;
	}
	os << endl;

	os << "Profile:" << endl;
	os << endl;

	for (ip = 0; ip < np_; ip++)
		os << "\t" << (acos(1.0 - (ip + 0.5) * pstep_) * 180 / PI);
	os << endl;

	for (ip = 0; ip < np_; ip++)
	{
		double f = 0;
		for (im = 0; im < nm_; im++)
			f += Dose(ip, im);
		os << "\t" << f / nm_;
	}
	os << endl;
}
