#include "mcScoreMatrixRZ.h"
#include "mcGeometry.h"
#include "mcTransport.h"
#include <float.h>

mcScoreMatrixRZ::mcScoreMatrixRZ(const char* module_name, int nThreads, int nr, int nz, double rmax, double zmin, double zmax)
	:mcScore(module_name, nThreads)
{
	m_nr = nr;
	m_nz = nz;
	m_rmax = rmax;
	m_rm2 = m_rmax * m_rmax;
	m_rstep = m_rmax / m_nr;
	m_zmin = zmin;
	m_zmax = zmax;
	m_zstep = (m_zmax - m_zmin) / m_nz;

	int len = nThreads * m_nr * m_nz;
	m_MAll = new double[len];
	memset(m_MAll, 0, len * sizeof(double));
	m_M = new double*[nThreads];
	m_M[0] = m_MAll;
	for (int i = 1; i < nThreads_; i++)
		m_M[i] = m_M[i - 1] + m_nr * m_nz;
}

mcScoreMatrixRZ::~mcScoreMatrixRZ()
{
	if (m_MAll) delete[] m_MAll;
	if (m_M) delete[] m_M;
}

void mcScoreMatrixRZ::ScoreFluence(const mcParticle& particle)
{ }

void mcScoreMatrixRZ::ScorePoint(double edep
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

void mcScoreMatrixRZ::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)
{
	double eddens = edep / (p1 - p0).length();

	// —ортируем точки так, чтобы первой была точка с меньшим Z.
	bool doChange = p0.z() > p1.z();
	geomVector3D vz1(doChange ? p1 : p0);
	geomVector3D vz2(doChange ? p0 : p1);
	geomVector3D v = vz2 - vz1;
	v.normalize();

	// ¬ектор направлени€, длина которого по Z равна 1.
	bool isOrto = fabs(v.z()) <= DBL_EPSILON;
	geomVector3D uz = isOrto ? geomVector3D() : v * (1. / v.z());

	// ќбрезаем трек плоскост€ми торцов цилиндра.
	int iz1 = (vz1.z() < m_zmin) ? -1 : int((vz1.z() - m_zmin) / m_zstep);
	int iz2 = (vz2.z() < m_zmin) ? -1 : int((vz2.z() - m_zmin) / m_zstep);

	// Scoring может находитьс€ внутри большого объекта и может потребоватьс€ 
	// отсекать частицы в принципе не попадающие в область интереса.
	if (iz2 < 0 || iz1 >= m_nz)
		return;

	if (iz1 < 0)
		vz1 += (vz2 - vz1) * ((m_zmin - vz1.z()) / (vz2.z() - vz1.z()) + DBL_EPSILON);
	if (iz2 >= m_nz)
		vz2 += (vz1 - vz2) * ((vz2.z() - m_zmax) / (vz1.z() - vz2.z()) + DBL_EPSILON);

	// ѕоскольку расчет пересечени€ с цилиндром задача более трудоемка€
	// и на практике пересечений цилиндров меньше пересечений плоскостей Z
	// (так как частицы преимущественно движутс€ вдоль Z),
	// то сначала определ€ем пересечени€ по радиусу.

	// »нтересно! Ѕыл прецендент, когда частица находилась в бесконечном слое так далеко,
	// что произошло переполнение числа int.
	// ѕоэтому положение по радиусу отсекаеи не по индексам, а по положению в пространстве.
	double rr = vz1.lengthXY();
	if (rr >= m_rmax) { if (vz2.lengthXY() >= m_rmax) return; }
	int ir = int(rr / m_rstep);
	if (ir > m_nr) ir = m_nr;

	while ((vz2 - vz1) * v > DBL_EPSILON)
	{
		double r = ir * m_rstep;
		double d;
		bool toCenter;

		// ќпредел€ем направление движени€ к центру или от него
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

		//  онец трека в цилиндрическом кольце
		geomVector3D vz12 = vz1 + (v * d);

		// ѕровер€ем не прошли ли мы дальнюю точку глобального трека.
		if ((vz12 - vz2) * v > 0)
			vz12 = vz2;

		// ѕеребираем слои по Z
		geomVector3D pz1 = vz1;
		int iz = int((pz1.z() - m_zmin) / m_zstep);

		while ((vz12 - pz1) * v > DBL_EPSILON)
		{
			if (isOrto) {
				ir = int(((vz12 + vz1)*0.5).lengthXY() / m_rstep);
				scoreEnergyInVoxel(iThread, ir, iz, eddens * (vz12 - vz1).length());
				break;
			}

			geomVector3D pz2 = vz1 + (uz * (m_zmin + m_zstep * (iz + 1) - vz1.z()));
			if ((pz2 - vz12) * v > 0)
				pz2 = vz12;

			// ¬ычисл€ем индекс кольца более надежным способом.
			// »з-за проблемы округлени€ на границе в противном случае возникают сбои индексов.
			ir = int(((pz2 + pz1)*0.5).lengthXY() / m_rstep);
			scoreEnergyInVoxel(iThread, ir, iz, eddens * (pz2 - pz1).length());

			pz1 = pz2;
			iz++;
		}

		vz1 = vz12;
		if (toCenter) ir--; else ir++;
	}
}

void mcScoreMatrixRZ::scoreEnergyInVoxel(int iThread, int ir, int iz, double edep)
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreMatrixRZ::scoreEnergyInVoxel: thread exceed container size");

	if (ir < m_nr && iz >= 0 && iz < m_nz) {
		m_M[iThread][ir*m_nz + iz] += edep;
		etotal_[iThread] += edep;
	}
}

void mcScoreMatrixRZ::CE2D()
{
	if (dconverted_)
		throw std::exception("mcScoreMatrixRZ::CE2D: Already converted to dose");
	dconverted_ = true;
	for (int iz = 0; iz < m_nz; iz++)
	{
		for (int ir = 0; ir < m_nr; ir++)
		{
			double r1 = ir*m_rstep;
			double r2 = r1 + m_rstep;
			double vol = PI * (r2*r2 - r1*r1)*m_zstep;
			for (int i = 0; i < nThreads_; i++)
				m_M[i][ir*m_nz + iz] /= vol;
		}
	}
}

double mcScoreMatrixRZ::Dose(int iThread, int ir, int iz) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreMatrixRZ::Dose: thread exceed container size");
	else if (iz < 0 || iz >= m_nz || ir < 0 || ir >= m_nr)
		return 0;
	return m_M[iThread][ir*m_nz + iz];
}

double mcScoreMatrixRZ::Dose(int ir, int iz) const
{
	double f = 0;
	if (iz >= 0 && iz < m_nz && ir >= 0 && ir < m_nr)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_M[i][ir*m_nz + iz];
	}
	return f;
}

void mcScoreMatrixRZ::dumpVRML(ostream& os) const
{
	os << "# mcScoreMatrixRZ Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int ir, iz, it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;

	os << "# mcScoreMatrixRZ Score: " << name_ << endl;
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
		//  онцентрические круги
		for (ir = 1; ir <= m_nr; ir++) {
			double r = m_rstep * ir;
			for (it = 0; it < 360; it += da) {
				geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), z) * mttow;
				os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
				p = geomVector3D(r*sin(mPi*(it + da)), r*cos(mPi*(it + da)), z) * mttow;
				os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
				count++;
			}
		}
		//// –адиусы
		//for(it=0; it<360; it+=da) {
		//  geomVector3D p = geomVector3D(0, 0, z) * mttow;
		//  os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
		//  p = geomVector3D(m_rmax*sin(mPi*it), m_rmax*cos(mPi*it), z) * mttow;
		//  os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
		//  count++;
		//}
	}

	//// Ћинии вдоль оси цилиндра
	//for(ir=1; ir<=m_nr; ir++) {
	//  double r = m_rstep * ir;
	//  for(it=0; it<360; it+=da) {
	//    geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), m_zmin) * mttow;
	//    os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	//    p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), m_zmax) * mttow;
	//    os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	//    count++;
	//  }
	//}

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	for (it = 0; it < count; it++)
		os << "      " << 2 * it << ' ' << 2 * it + 1 << " -1" << endl;
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;
}

void mcScoreMatrixRZ::dumpStatistic(ostream& os) const
{
	int ir, iz;
	int nr = Nr();
	int nz = Nz();
	double rmax = MaxR();
	double zmin = MinZ();
	double zmax = MaxZ();
	double rstep = StepR();
	double zstep = StepZ();

	mcScore::dumpStatistic(os);

	if (dconverted_)
		os << "Dose matrix" << endl;
	else
		os << "Energy deposition matrix" << endl;

	os << "NR\tNz\tRMAX\tZ1\tZ2" << endl;
	os << nr << '\t' << nz << '\t' << rmax << '\t' << zmin << '\t' << zmax << endl << endl;

	//for( iz=0; iz < nz; iz++ )
	//	os << '\t' << zmin + (double(iz) + 0.5) * zstep;
	//os << endl;

	//for( ir=0; ir < nr; ir++ )
	//{
	//	os << (double(ir) + 0.5) * rstep;

	//    for( iz=0; iz < nz; iz++ )
	//	    os << '\t' << m[ir*nz + iz];
	//    os << endl;
	//}

	// this code - for save transpose matrix
	for (ir = 0; ir < nr; ir++)
		os << '\t' << (double(ir) + 0.5) * rstep;
	os << endl;

	for (iz = 0; iz < nz; iz++)
	{
		os << zmin + (double(iz) + 0.5) * zstep;
		for (ir = 0; ir < nr; ir++)
			os << '\t' << Dose(ir, iz);
		os << endl;
	}
}
