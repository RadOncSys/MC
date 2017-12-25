#include "mcScoreMatrixCylinderAzimut.h"
#include "mcGeometry.h"
#include "mcTransport.h"
#include <float.h>

mcScoreMatrixCylinderAzimut::mcScoreMatrixCylinderAzimut(const char* module_name, int nThreads, int nr, int na, double rmax, double zmin, double zmax)
	:mcScore(module_name, nThreads)
{
	m_nr = nr;
	m_na = na;
	m_rmax = rmax;
	m_rm2 = m_rmax * m_rmax;
	m_r2step = m_rm2 / m_nr;
	m_zmin = zmin;
	m_zmax = zmax;
	m_astep = 2 * PI / m_na;

	int len = nThreads * m_nr * m_na;
	m_MAll = new double[len];
	memset(m_MAll, 0, len * sizeof(double));
	m_M = new double*[nThreads];
	m_M[0] = m_MAll;
	for (int i = 1; i < nThreads_; i++)
		m_M[i] = m_M[i - 1] + m_nr * m_na;
}

mcScoreMatrixCylinderAzimut::~mcScoreMatrixCylinderAzimut()
{
	if (m_MAll) delete[] m_MAll;
	if (m_M) delete[] m_M;
}

void mcScoreMatrixCylinderAzimut::ScoreFluence(const mcParticle& particle)
{ }

void mcScoreMatrixCylinderAzimut::ScorePoint(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p)
{
	if (p.z() < m_zmin || p.z() > m_zmax)
		return;

	int ir, ia;
	getVoxelAtPoint(p, ir, ia);
	if (ir < 0)
		return;

	scoreEnergyInVoxel(iThread, ir, ia, edep);
}

void mcScoreMatrixCylinderAzimut::ScoreLine(double edep
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

	// “рек за пределами торцевых плоскостей скоринга целиком
	if (vz2.z() <= m_zmin || vz1.z() >= m_zmax)
		return;

	// ќбрезаем трек плоскост€ми торцов цилиндра.
	if (vz1.z() < m_zmin)
		vz1 += v * ((m_zmin - vz1.z()) / v.z());
	if (vz2.z() > m_zmax)
		vz2 += v * ((m_zmax - vz2.z()) / v.z());

	// — этого момента задача становитс€ плоской.
	// Ќужно определить длины треков в каждой €чейке, определ€емой кругами и радиусами.

	// ѕоворачиваем направление движени€ против часовой стрелке (увеличение угла).
	doChange = (vz2 ^ vz1).z() > 0;
	if (doChange)
	{
		geomVector3D vtmp(vz1);
		vz1 = vz2;
		vz2 = vtmp;
		v *= -1;
	}

	// –ассчитываем длины треков в €чейках просто двига€сь вдоль траетории
	// и определ€€ точки пересечени€ с границами €чеек.
	while (true)
	{
		int ir, ia;
		getVoxelAtPoint(vz1, ir, ia);
		bool isMoveInside = (vz1.x() * v.x() + vz1.y() * v.y()) < 0;
		if (ir < 0 && !isMoveInside) break;

		// ≈сли снаружи, то сначала определ€ем точку входа в скоринг.
		if (ir < 0)
		{
			double d = mcGeometry::getDistanceToInfiniteCylinderOutside(vz1, v, m_rmax);
			if (d == DBL_MAX) break;	// вообще пролетели мимо скоринга
			d += FLT_EPSILON;			// гарантируем пересечение границы
			if (d >= (vz2 - vz1).length()) break;
			vz1 += v * d;
			continue;					// обрезали внешнюю часть трека и повтор€ем движение уже реально по вокселам 
		}
		
		// “еперь вопрос что пересечем раньше, окружность или радиус.
		double dr = isMoveInside ? 
			mcGeometry::getDistanceToInfiniteCylinderOutside(vz1, v, sqrt(ir*m_r2step)) :
			mcGeometry::getDistanceToInfiniteCylinderInside(vz1, v, sqrt((ir + 1)*m_r2step));
		
		// Ќормаль к радиальной плоскости
		double a = (ia + 1) * m_astep;
		geomVector3D n(-sin(a), cos(a), 0);
		double h = -(vz1 * n);
		double f = v * n;	// косинус между нормалью и направление трека
		double da = f < FLT_EPSILON ? DBL_MAX : h / f;	// отрицательное f означает вылетание из конуса без пересечени€ его границ
		double d = MIN(da, dr) + FLT_EPSILON;
		double dmax = (vz2 - vz1).length();

		if (d >= dmax)
		{
			// “рек заканчиваетс€ в данном вокселе
			scoreEnergyInVoxel(iThread, ir, ia, dmax * eddens);
			break;
		}
		scoreEnergyInVoxel(iThread, ir, ia, d * eddens);
		vz1 += v * d;
	}
}

void mcScoreMatrixCylinderAzimut::scoreEnergyInVoxel(int iThread, int ir, int ia, double edep)
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreMatrixCylinderAzimut::scoreEnergyInVoxel: thread exceed container size");

	if (ir < m_nr && ia >= 0 && ia < m_na) {
		m_M[iThread][ir*m_na + ia] += edep;
		etotal_[iThread] += edep;
	}
}

void mcScoreMatrixCylinderAzimut::CE2D()
{
	if (dconverted_)
		throw std::exception("mcScoreMatrixCylinderAzimut::CE2D: Already converted to dose");
	dconverted_ = true;
	for (int ia = 0; ia < m_na; ia++)
	{
		for (int ir = 0; ir < m_nr; ir++)
		{
			double vol = PI * m_r2step / m_na;
			for (int i = 0; i < nThreads_; i++)
				m_M[i][ir*m_na + ia] /= vol;
		}
	}
}

double mcScoreMatrixCylinderAzimut::Dose(int iThread, int ir, int ia) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreMatrixCylinderAzimut::Dose: thread exceed container size");
	else if (ia < 0 || ia >= m_na || ir < 0 || ir >= m_nr)
		return 0;
	return m_M[iThread][ir*m_na + ia];
}

double mcScoreMatrixCylinderAzimut::Dose(int ir, int ia) const
{
	double f = 0;
	if (ia >= 0 && ia < m_na && ir >= 0 && ir < m_nr)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_M[i][ir*m_na + ia];
	}
	return f;
}

void mcScoreMatrixCylinderAzimut::getVoxelAtPoint(const geomVector3D& p, int& ir, int& ia)
{
	ir = -1; ia = -1;
	double r2 = p.sqLengthXY();
	if (r2 >= m_rm2) return;
	double r = sqrt(r2);

	double x = p.x(), y = p.y();
	double angle = acos(x / r);
	if (y < 0) angle = 2 * PI - angle;
	ia = (int)(angle / m_astep);
	if (ia >= m_na) {
		ia = -1;
		return;
	}

	ir = int(m_nr * r2 / m_rm2);
	if (ir >= m_nr) {
		ir = -1; 
		return;
	}
}

void mcScoreMatrixCylinderAzimut::dumpVRML(ostream& os) const
{
	os << "# mcScoreMatrixCylinderAzimut Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int ir, ia, it, count = 0;
	double mPi = PI / 180;

	os << "# mcScoreMatrixCylinderAzimut Score: " << name_ << endl;
	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	for (int iz = 0; iz < 2; iz++)
	{
		double z = iz == 0 ? m_zmin : m_zmax;

		//  онцентрические круги
		for (ir = 1; ir <= m_nr; ir++) 
		{
			double r = sqrt(m_r2step * ir);
			for (ia = 0; ia <= m_na; ia ++)
			{
				geomVector3D p = geomVector3D(r*sin(m_astep*ia), r*cos(m_astep*ia), z) * mttow;
				os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
				p = geomVector3D(r*sin(m_astep*(ia+1)), r*cos(m_astep*(ia + 1)), z) * mttow;
				os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
				count++;
			}
		}
		// –адиусы
		for (ia = 0; ia < m_na; ia++)
		{
			geomVector3D p = geomVector3D(0, 0, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(m_rmax*sin(m_astep*ia), m_rmax*cos(m_astep*ia), z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			count++;
		}
	}

	// Ћинии вдоль оси цилиндра
	for(ir=1; ir<=m_nr; ir++)
	{
		double r = sqrt(m_r2step * ir);
		for (ia = 0; ia < m_na; ia++)
		{
			geomVector3D p = geomVector3D(r*sin(m_astep*ia), r*cos(m_astep*ia), m_zmin) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r*sin(m_astep*ia), r*cos(m_astep*ia), m_zmax) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			count++;
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

void mcScoreMatrixCylinderAzimut::dumpStatistic(ostream& os) const
{
	int ir, ia;
	int nr = Nr();
	int na = Na();
	double rmax = MaxR();
	double zmin = MinZ();
	double zmax = MaxZ();
	double r2step = StepR2();
	double astep = StepA();

	mcScore::dumpStatistic(os);

	if (dconverted_)
		os << "Dose matrix" << endl;
	else
		os << "Energy deposition matrix" << endl;

	os << "NR\tNa\tRMAX\tZ1\tZ2" << endl;
	os << nr << '\t' << na << '\t' << rmax << '\t' << zmin << '\t' << zmax << endl << endl;

	// this code - for save transpose matrix
	for (ir = 0; ir < nr; ir++)
		os << '\t' << sqrt((ir + 0.5) * r2step);
	os << endl;

	for (ia = 0; ia < na; ia++)
	{
		os << (ia + 0.5) * astep;
		for (ir = 0; ir < nr; ir++)
			os << '\t' << Dose(ir, ia);
		os << endl;
	}
}
