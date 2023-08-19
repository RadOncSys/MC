#include "mcScoreBrachy.h"
#include "mcGeometry.h"
#include "mcTransport.h"
#include "ProfileProcessor.h"
#include <float.h>

mcScoreBrachy::mcScoreBrachy(const char* module_name, int nThreads, double sourceLength)
	:mcScore(module_name, nThreads), L_(sourceLength)
{
	// Учитывая что в статистике может применяться сглаживание не стоит важную область делать пограничной.
	// Поэтому в ее сборе частично заходим в область отрицательных Z.
	// Сетку по Z подбираем так, чтобы не нужно было интерполировать.

	m_nr = 320;
	m_nz = 441;
	m_rmax = 16;
	m_rm2 = m_rmax * m_rmax;
	m_rstep = m_rmax / m_nr;
	m_zmin = -11.025;
	m_zmax = 11.025;
	m_zstep = (m_zmax - m_zmin) / m_nz;

	int len = nThreads * m_nr * m_nz;
	m_MAll = new double[len];
	memset(m_MAll, 0, len * sizeof(double));
	m_M = new double* [nThreads];
	m_M[0] = m_MAll;
	for (int i = 1; i < nThreads_; i++)
		m_M[i] = m_M[i - 1] + m_nr * m_nz;
}

mcScoreBrachy::~mcScoreBrachy()
{
	if (m_MAll) delete[] m_MAll;
	if (m_M) delete[] m_M;
}

void mcScoreBrachy::ScoreFluence(const mcParticle& particle)
{ }

void mcScoreBrachy::ScorePoint(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p)
{
	if (p.z() < m_zmin)
		return;
	int iz = int((p.z() - m_zmin) / m_zstep);
	if (iz >= m_nz)
		return;

	double r = p.lengthXY();
	int ir = int(r / m_rstep);
	if (ir >= m_nr)
		return;

	scoreEnergyInVoxel(iThread, ir, iz, edep);
}

void mcScoreBrachy::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)
{
	double eddens = edep / (p1 - p0).length();

	// Сортируем точки так, чтобы первой была точка с меньшим Z.
	bool doChange = p0.z() > p1.z();
	geomVector3D vz1(doChange ? p1 : p0);
	geomVector3D vz2(doChange ? p0 : p1);
	geomVector3D v = vz2 - vz1;
	v.normalize();

	// Вектор направления, длина которого по Z равна 1.
	bool isOrto = fabs(v.z()) <= DBL_EPSILON;
	geomVector3D uz = isOrto ? geomVector3D() : v * (1. / v.z());

	// Обрезаем трек плоскостями торцов цилиндра.
	int iz1 = (vz1.z() < m_zmin) ? -1 : int((vz1.z() - m_zmin) / m_zstep);
	int iz2 = (vz2.z() < m_zmin) ? -1 : int((vz2.z() - m_zmin) / m_zstep);

	// Scoring может находиться внутри большого объекта и может потребоваться 
	// отсекать частицы в принципе не попадающие в область интереса.
	if (iz2 < 0 || iz1 >= m_nz)
		return;

	if (iz1 < 0)
		vz1 += (vz2 - vz1) * ((m_zmin - vz1.z()) / (vz2.z() - vz1.z()) + DBL_EPSILON);
	if (iz2 >= m_nz)
		vz2 += (vz1 - vz2) * ((vz2.z() - m_zmax) / (vz1.z() - vz2.z()) + DBL_EPSILON);

	// Поскольку расчет пересечения с цилиндром задача более трудоемкая
	// и на практике пересечений цилиндров меньше пересечений плоскостей Z
	// (так как частицы преимущественно движутся вдоль Z),
	// то сначала определяем пересечения по радиусу.

	// Интересно! Был прецендент, когда частица находилась в бесконечном слое так далеко,
	// что произошло переполнение числа int.
	// Поэтому положение по радиусу отсекаеи не по индексам, а по положению в пространстве.
	double rr = vz1.lengthXY();
	if (rr >= m_rmax) { if (vz2.lengthXY() >= m_rmax) return; }
	int ir = int(rr / m_rstep);
	if (ir > m_nr) ir = m_nr;

	while ((vz2 - vz1) * v > DBL_EPSILON)
	{
		double r = ir * m_rstep;
		double d;
		bool toCenter;

		// Определяем направление движения к центру или от него
		if (ir > 0 && vz1 * v < 0) // к центру
		{
			d = mcGeometry::getDistanceToInfiniteCylinderOutside(vz1, v, r);
			toCenter = true;
			if (d == DBL_MAX) {    // мимо внутреннего цилиндра
				if (ir >= m_nr) return;
				d = mcGeometry::getDistanceToInfiniteCylinderInside(vz1, v, r + m_rstep);
				toCenter = false;
			}
			else {
				if (ir >= m_nr) {
					vz1 += v * d;
					ir--;
					continue;
				}
			}
		}
		else
		{
			if (ir >= m_nr) return;
			d = mcGeometry::getDistanceToInfiniteCylinderInside(vz1, v, r + m_rstep);
			toCenter = false;
		}

		// Конец трека в цилиндрическом кольце
		geomVector3D vz12 = vz1 + (v * d);

		// Проверяем не прошли ли мы дальнюю точку глобального трека.
		if ((vz12 - vz2) * v > 0)
			vz12 = vz2;

		// Перебираем слои по Z
		geomVector3D pz1 = vz1;
		int iz = int((pz1.z() - m_zmin) / m_zstep);

		while ((vz12 - pz1) * v > DBL_EPSILON)
		{
			if (isOrto) {
				ir = int(((vz12 + vz1) * 0.5).lengthXY() / m_rstep);
				scoreEnergyInVoxel(iThread, ir, iz, eddens * (vz12 - vz1).length());
				break;
			}

			geomVector3D pz2 = vz1 + (uz * (m_zmin + m_zstep * (iz + 1) - vz1.z()));
			if ((pz2 - vz12) * v > 0)
				pz2 = vz12;

			// Вычисляем индекс кольца более надежным способом.
			// Из-за проблемы округления на границе в противном случае возникают сбои индексов.
			ir = int(((pz2 + pz1) * 0.5).lengthXY() / m_rstep);
			scoreEnergyInVoxel(iThread, ir, iz, eddens * (pz2 - pz1).length());

			pz1 = pz2;
			iz++;
		}

		vz1 = vz12;
		if (toCenter) ir--; else ir++;
	}
}

void mcScoreBrachy::scoreEnergyInVoxel(int iThread, int ir, int iz, double edep)
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBrachy::scoreEnergyInVoxel: thread exceed container size");

	if (ir < m_nr && iz >= 0 && iz < m_nz) {
		m_M[iThread][ir * m_nz + iz] += edep;
		etotal_[iThread] += edep;
	}
}

void mcScoreBrachy::CE2D()
{
	if (dconverted_)
		throw std::exception("mcScoreBrachy::CE2D: Already converted to dose");
	dconverted_ = true;
	for (int iz = 0; iz < m_nz; iz++)
	{
		for (int ir = 0; ir < m_nr; ir++)
		{
			double r1 = ir * m_rstep;
			double r2 = r1 + m_rstep;
			double vol = PI * (r2 * r2 - r1 * r1) * m_zstep;
			for (int i = 0; i < nThreads_; i++)
				m_M[i][ir * m_nz + iz] /= vol;
		}
	}
}

double mcScoreBrachy::Dose(int iThread, int ir, int iz) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBrachy::Dose: thread exceed container size");
	else if (iz < 0 || iz >= m_nz || ir < 0 || ir >= m_nr)
		return 0;
	return m_M[iThread][ir * m_nz + iz];
}

double mcScoreBrachy::Dose(int ir, int iz) const
{
	double f = 0;
	if (iz >= 0 && iz < m_nz && ir >= 0 && ir < m_nr)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_M[i][ir * m_nz + iz];
	}
	return f;
}

void mcScoreBrachy::dumpVRML(ostream& os) const
{
	os << "# mcScoreBrachy Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int ir, iz, it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;

	os << "# mcScoreBrachy Score: " << name_ << endl;
	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	for (iz = 0; iz <= m_nz; iz++)
	{
		double z = m_zmin + iz * m_zstep;
		// Концентрические круги
		for (ir = 1; ir <= m_nr; ir++) 
		{
			// Рисуем только линии на поверхности цилиндра детектирующей матрицы
			if (iz == 0 || iz == m_nz || ir == m_nr)
			{
				double r = m_rstep * ir;
				for (it = 0; it < 360; it += da) {
					geomVector3D p = geomVector3D(r * sin(mPi * it), r * cos(mPi * it), z) * mttow;
					os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
					p = geomVector3D(r * sin(mPi * (it + da)), r * cos(mPi * (it + da)), z) * mttow;
					os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
					count++;
				}
			}
		}
	}

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	for (it = 0; it < count; it++)
		os << "      " << 2 * it << ' ' << 2 * it + 1 << " -1" << endl;
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;
}

void mcScoreBrachy::dumpStatistic(ostream& os) const
{
	int ir, iz;
	mcScore::dumpStatistic(os);

	if (!dconverted_)
		throw std::exception("mcScoreBrachy::dumpStatistic: Energy deposition must be converted to Dose before dump statistic");

	// Создаем массив для установки матрицы F(r,thetta)

	auto D = std::make_unique<std::vector<std::vector<double>>>(m_nr, std::vector<double>(m_nz, 0));
	auto& dref = *D.get();
	for (ir = 0; ir < m_nr; ir++)
		for (iz = 0; iz < m_nz; iz++)
			dref[ir][iz] = Dose(ir, iz);

	//auto F = ProfileProcessor::SmoothSG2D(dref);
	//auto& fref = *F.get();

	// g(r)
	double r_t[] = { 0.1, 0.2, 0.3, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
	int nr = sizeof(r_t) / sizeof(double);

	os << "Radial dose function, g(r)" << endl;
	os << "--------------------------" << endl;
	os << "r[cm]\tg(r)" << endl;

	int idz0 = int(-m_zmin / m_zstep);
	int idr = int(1.0 / m_rstep - 0.5);
	double f = 1.0 / m_rstep - idr - 0.5;
	//double gr0 = fref[idr][idz0] * (1 - f) + fref[idr + 1][idz0] * f;
	double gr0 = dref[idr][idz0] * (1 - f) + dref[idr + 1][idz0] * f;
	double G0 = 2 * atan(L_ / 2) / L_;

	for (ir = 0; ir < nr; ir++)
	{
		// По радиусу придется интерполировать, так как он начинается с 0
		// и нужные точки оказываются между пикселами.
		double r = r_t[ir];
		idr = int(r / m_rstep - 0.5);
		f = r / m_rstep - idr - 0.5;
		double G = 2 * atan(L_ / (2 * r)) / (L_ * r);
		//double gr = ((fref[idr][idz0] * (1 - f) + fref[idr + 1][idz0] * f) / gr0) * (G0 / G);
		//os << r << "\t" << gr;
		double gr = ((dref[idr][idz0] * (1 - f) + dref[idr + 1][idz0] * f) / gr0) * (G0 / G);
		os << r << "\t" << gr << endl;
	}
	os << std::endl;

	// Шаблоны самплинга двумерной функции анизотропии F(r, theta)
	double fr_t[] = { 0.25, 0.5, 1, 2, 3, 5 };
	double fa_t[] = { 0, 1, 2, 3, 5, 7, 9, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 75, 90, 105, 120, 130, 135, 140, 145, 150, 155, 160, 165, 168, 171, 173, 175	};

	int nfrt = sizeof(fr_t) / sizeof(double);
	int nfat = sizeof(fa_t) / sizeof(double);

	os << "2D anisotropy function, F(r, thetta)" << endl;
	os << "------------------------------------" << endl;
	os << std::endl;

	for (ir = 0; ir < nfrt; ir++)
		os << '\t' << fr_t[ir];
	os << std::endl;

	for (iz = 0; iz < nfat; iz++)
	{
		double a = PI * fa_t[iz] / 180;
		os << fa_t[iz];

		for (ir = 0; ir < nfrt; ir++)
		{
			double r = fr_t[ir];
			double x = r * cos(a);
			double y = r * sin(a);

			double G = (fa_t[iz] == 0) ? L_ / (r * r - L_ * L_ / 4) :
				(atan((x - L_ / 2) / y) - atan((x + L_ / 2) / y)) / y;
			if (G < 0) G = -G;
			double G0 = 2 * atan(L_ / (2 * r)) / r;

			// D0
			int idy0 = int(r / m_rstep - 0.5);
			f = r / m_rstep - idy0 - 0.5;
			double d0 = dref[idy0][idz0] * (1 - f) + dref[idy0 + 1][idz0] * f;

			// D
			int idx = int((x - m_zmin) / m_zstep - 0.5);
			int idy = int(y / m_rstep - 0.5);

			double d00 = dref[idy][idx];
			double d01 = dref[idy][idx + 1];
			double d10 = dref[idy + 1][idx];
			double d11 = dref[idy + 1][idx + 1];

			double fx = (x - m_zmin) / m_zstep - idx - 0.5;
			double fy = idy == 0 ? 0 : y / m_rstep - idy - 0.5;

			double d = (d00 * (1 - fx) + d01 * fx) * (1 - fy) + (d10 * (1 - fx) + d11 * fx) * fy;

			double frt = (d * G0) / (d0 * G);

			os << "\t" << frt;
		}
		os << std::endl;
	}
	os << std::endl;

	//// Dose matrix
	//os << std::endl;
	//for (ir = 0; ir < m_nr; ir += 5)
	//	os << '\t' << (double(ir) + 0.5) * m_rstep;
	//os << endl;

	//for (iz = 0; iz < m_nz; iz += 5)
	//{
	//	os << m_zmin + (double(iz) + 0.5) * m_zstep;
	//	for (ir = 0; ir < m_nr; ir += 5)
	//		os << '\t' << dref[ir][iz];
	//	os << endl;
	//}
}
