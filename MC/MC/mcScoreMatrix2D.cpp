#include "mcScoreMatrix2D.h"
#include "mcGeometry.h"
#include "mcTransport.h"

mcScoreMatrix2D::mcScoreMatrix2D(const char* module_name, int nThreads, int nx, int nz, const geomVector3D& v1, const geomVector3D& v2)
	:mcScore(module_name, nThreads)
	, m_nx(nx), m_nz(nz)
{
	if (v1.y() > v2.y())
	{
		m_ymin = v2.y();
		m_ymax = v1.y();
	}
	else
	{
		m_ymin = v1.y();
		m_ymax = v2.y();
	}
	m_xmin = v1.x();
	m_xmax = v2.x();

	m_zmin = v1.z();
	m_zmax = v2.z();

	m_xstep = (m_xmax - m_xmin) / m_nx;
	m_zstep = (m_zmax - m_zmin) / m_nz;

	int len = nThreads * m_nx * m_nz;
	m_MAll = new double[len];
	memset(m_MAll, 0, len * sizeof(double));
	m_M = new double*[nThreads];
	m_M[0] = m_MAll;
	for (int i = 1; i < nThreads_; i++)
		m_M[i] = m_M[i - 1] + m_nx * m_nz;

}

mcScoreMatrix2D::~mcScoreMatrix2D()
{
	if (m_MAll) delete[] m_MAll;
	if (m_M) delete[] m_M;
}

void mcScoreMatrix2D::ScoreFluence(const mcParticle& particle)
{ }

void mcScoreMatrix2D::ScorePoint(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p)
{
	if (p.y() < m_ymin || p.y() >= m_ymax) return;
	if (p.x() <= m_xmin || p.z() <= m_zmin)
		return;

	int ix = int((p.x() - m_xmin) / m_xstep);
	int iz = int((p.z() - m_zmin) / m_zstep);
	if (ix >= m_nx || iz >= m_nz)
		return;

	m_M[iThread][iz*m_nx + ix] += edep;
	etotal_[iThread] += edep;
}

void mcScoreMatrix2D::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)
{
	double x1 = p0.x(), y1 = p0.y(), z1 = p0.z();
	double x2 = p1.x(), y2 = p1.y(), z2 = p1.z();

	// Make program faster
	if (x1 < m_xmin && x2 < m_xmin) return;
	if (x1 > m_xmax && x2 > m_xmax) return;

	if (y1 < m_ymin && y2 < m_ymin) return;
	if (y1 > m_ymax && y2 > m_ymax) return;

	if (z1 < m_zmin && z2 < m_zmin) return;
	if (z1 > m_zmax && z2 > m_zmax) return;

	// Index limits
	int ix1 = (x1 < m_xmin) ? -1 : int((x1 - m_xmin) / m_xstep);
	int ix2 = (x2 < m_xmin) ? -1 : int((x2 - m_xmin) / m_xstep);
	int iz1 = (z1 < m_zmin) ? -1 : int((z1 - m_zmin) / m_zstep);
	int iz2 = (z2 < m_zmin) ? -1 : int((z2 - m_zmin) / m_zstep);

	// When path is completely within one voxel
	if (ix1 == ix2 && iz1 == iz2 &&
		y1 >= m_ymin && y1 <= m_ymax && y2 >= m_ymin && y2 <= m_ymax)
	{
		m_M[iThread][iz1*m_nx + ix1] += edep;
		etotal_[iThread] += edep;
		return;
	}

	ix1 = (ix1 < 0) ? 0 : (ix1 >= m_nx) ? m_nx - 1 : ix1;
	ix2 = (ix2 < 0) ? 0 : (ix2 >= m_nx) ? m_nx - 1 : ix2;
	iz1 = (iz1 < 0) ? 0 : (iz1 >= m_nz) ? m_nz - 1 : iz1;
	iz2 = (iz2 < 0) ? 0 : (iz2 >= m_nz) ? m_nz - 1 : iz2;

	int i;
	if (ix1 > ix2) { i = ix1; ix1 = ix2; ix2 = i; }
	if (iz1 > iz2) { i = iz1; iz1 = iz2; iz2 = i; }

	// Energy lose per unit path.
	double f = edep / (p1 - p0).length();

	// Go trough voxels
	double yv1 = m_ymin, yv2 = m_ymax;
	for (int iz = iz1; iz <= iz2; iz++)
	{
		double zv1 = m_zmin + m_zstep * iz, zv2 = zv1 + m_zstep;
		for (int ix = ix1; ix <= ix2; ix++)
		{
			double xv1 = m_xmin + m_xstep * ix, xv2 = xv1 + m_xstep;
			double trl = mcGeometry::TrLenInVoxel(x1, y1, z1, x2, y2, z2, xv1, yv1, zv1, xv2, yv2, zv2);
			m_M[iThread][iz*m_nx + ix] += f * trl;
			etotal_[iThread] += f * trl;
		}
	}
	return;
}

void mcScoreMatrix2D::CE2D()
{
	if (dconverted_)
		throw std::exception("mcScoreMatrix2D::CE2D: Already converted to dose");
	dconverted_ = true;
	double f = m_nx * m_nz / (density_ * (m_xmax - m_xmin) * (m_ymax - m_ymin) * (m_zmax - m_zmin));
	int ix, iz;
	for (iz = 0; iz < m_nz; iz++) {
		for (ix = 0; ix < m_nx; ix++)
			for (int i = 0; i < nThreads_; i++)
				m_M[i][iz*m_nx + ix] *= f;
	}
}

void mcScoreMatrix2D::mul(double f)
{
	int ix, iz;
	for (iz = 0; iz < m_nz; iz++) {
		for (ix = 0; ix < m_nx; ix++)
			for (int i = 0; i < nThreads_; i++)
				m_M[i][iz*m_nx + ix] *= f;
	}
}

double mcScoreMatrix2D::Dose(int iThread, int ix, int iz) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreMatrix2D::Dose: thread exceed container size");
	else if (iz < 0 || iz >= m_nz || ix < 0 || ix >= m_nx)
		return 0;
	return m_M[iThread][iz*m_nx + ix];
}

double mcScoreMatrix2D::Dose(int ix, int iz) const
{
	double f = 0;
	if (iz >= 0 && iz < m_nz && ix >= 0 && ix < m_nx)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_M[i][iz*m_nx + ix];
	}
	return f;
}


double mcScoreMatrix2D::Sigma(int ix, int iz) const
{
	double f = Dose(ix, iz);// суммарна€ доза по батчам (средн€€ доза умнож на кол-во потоков)
	if (f == 0)
		return 0;
	f /= nThreads_; // средн€€ доза по батчам
	double sigma = 0;
	for (int register i = 0; i < nThreads_; i++)///цикл дл€ суммировани€ сигма по батчам
	{
		sigma += pow((m_M[i][iz*m_nx + ix] - f), 2);
	}
	sigma = sqrt(sigma / nThreads_);
	return sigma / f;
}


void mcScoreMatrix2D::dumpStatistic(ostream& os) const
{
	int ix, iz;
	int nx = Nx();
	int nz = Nz();
	double xmin = MinX();
	double xmax = MaxX();
	double ymin = MinY();
	double ymax = MaxY();
	double zmin = MinZ();
	double zmax = MaxZ();
	double xstep = StepX();
	double zstep = StepZ();

	mcScore::dumpStatistic(os);

	if (dconverted_)
		os << "Energy:" << endl;
	else
		os << "Dose:" << endl;

	os << "NX\tNZ\tX1\tY1\tZ1\tX2\tY2\tZ2" << endl;
	os << nx << '\t' << nz << '\t' << xmin << '\t' << ymin << '\t' << zmin << '\t'
		<< xmax << '\t' << ymax << '\t' << zmax << endl << endl;

	for (ix = 0; ix < nx; ix++)
		os << '\t' << xmin + (double(ix) + 0.5) * xstep;
	os << endl;

	for (iz = 0; iz < nz; iz++)
	{
		os << zmin + (double(iz) + 0.5) * zstep;
		for (ix = 0; ix < nx; ix++)
			os << '\t' << Dose(ix, iz);
		os << endl;
	}

	// this code - for Sigma
	os << endl;
	os << "Sigma matrix" << endl;
	os << endl;
	for (ix = 0; ix < nx; ix++)
		os << '\t' << (double(ix) + 0.5) * xstep;
	os << endl;

	for (iz = 0; iz < nz; iz++)
	{
		os << zmin + (double(iz) + 0.5) * zstep;
		for (ix = 0; ix < nx; ix++)
			os << '\t' << Sigma(ix, iz);
		os << endl;
	}

	double sigmaMean = 0;
	int count = 0;
	for (iz = 4 * nz / 10; iz < 6 * nz / 10; iz++)
	{
		for (ix = 3 * nx / 10; ix < 7 * nx / 10; ix++)
		{
			sigmaMean += Sigma(ix, iz);
			count++;
		}
	}
	os << endl;
	os << "Mean sigma = \t" << sigmaMean / count << endl;
	os << endl;

}



void mcScoreMatrix2D::dumpVRML(ostream& os)const
{
	os << "# Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	// TODO: Ќадо задействовать преобразование координат в мировую систему !!!

	unsigned i, j, k, ni = m_nx + 1, nj = 2, nk = m_nz + 1;

	os << "# Score: " << name_ << endl;
	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "      transparency " << transparancy_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	// √рани XY (четные точки нам Zmin, нечетные - Zmax)
	unsigned ncount = 0;
	for (j = 0; j < nj; j++) {
		double y = m_ymin + j*(m_ymax - m_ymin) / (nj - 1);
		for (i = 0; i < ni; i++) {
			double x = m_xmin + i*m_xstep;
			geomVector3D p = geomVector3D(x, y, m_zmin) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(x, y, m_zmax) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			ncount++;
		}
	}
	// √рани XZ
	for (k = 0; k < nk; k++) {
		double z = m_zmin + k*m_zstep;
		for (i = 0; i < ni; i++) {
			double x = m_xmin + i*m_xstep;
			geomVector3D p = geomVector3D(x, m_ymin, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(x, m_ymax, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			ncount++;
		}
	}
	// √рани YZ
	for (k = 0; k < nk; k++) {
		double z = m_zmin + k*m_zstep;
		for (j = 0; j < nj; j++) {
			double y = m_ymin + j*(m_ymax - m_ymin) / (nj - 1);
			geomVector3D p = geomVector3D(m_xmin, y, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(m_xmax, y, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			ncount++;
		}
	}

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	for (i = 0; i < ncount; i++)
		os << "      " << 2 * i << ' ' << 2 * i + 1 << " -1" << endl;
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;
}


