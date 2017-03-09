#include "mcScoreConicalRZ.h"
#include "mcGeometry.h"
#include "mcTransport.h"
#include <float.h>
#include "ProfileProcessor.h"
#include <strstream>

mcScoreConicalRZ::mcScoreConicalRZ(const char* module_name, int nThreads, int nr, int nz, double rmax,
	double zmin, double zmax, double ziso, double sad)
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
	m_iso = ziso;
	m_sad = sad;

	int len = nThreads * m_nr * m_nz;
	m_MAll = new double[len];
	memset(m_MAll, 0, len * sizeof(double));
	m_M = new double*[nThreads];
	m_M[0] = m_MAll;
	for (int i = 1; i < nThreads_; i++)
		m_M[i] = m_M[i - 1] + m_nr * m_nz;

	//Sigma = new double[nThreads];
}

mcScoreConicalRZ::~mcScoreConicalRZ()
{
	if (m_MAll) delete[] m_MAll;
	if (m_M) delete[] m_M;
}

void mcScoreConicalRZ::ScoreFluence(const mcParticle& particle)
{ }

void mcScoreConicalRZ::ScorePoint(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p)
{
	double z = p.z();
	if (z < m_zmin)
		return;
	int iz = int((z - m_zmin) / m_zstep);
	if (iz >= m_nz)
		return;

	double r = p.lengthXY() / (1.0 + (z - m_iso) / m_sad);
	int ir = int(r / m_rstep);
	if (ir >= m_nr)
		return;

	scoreEnergyInVoxel(iThread, ir, iz, edep);
}

void mcScoreConicalRZ::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)
{
	// ��������� ����� ���, ����� ������ ���� ����� � ������� Z.
	bool doChange = p0.z() > p1.z();
	geomVector3D vz1(doChange ? p1 : p0);
	geomVector3D vz2(doChange ? p0 : p1);
	geomVector3D v = vz2 - vz1;
	v.normalize();

	// �������� ���� ����������� ������ ��������.
	int iz1 = (vz1.z() < m_zmin) ? -1 : int((vz1.z() - m_zmin) / m_zstep);
	int iz2 = (vz2.z() < m_zmin) ? -1 : int((vz2.z() - m_zmin) / m_zstep);

	// Scoring ����� ���������� ������ �������� ������� � ����� ������������� 
	// �������� ������� � �������� �� ���������� � ������� ��������.
	if (iz2 < 0 || iz1 >= m_nz)
		return;

	if (iz1 < 0)
		vz1 += (vz2 - vz1) * ((m_zmin - vz1.z()) / (vz2.z() - vz1.z()) + DBL_EPSILON);
	if (iz2 >= m_nz)
		vz2 += (vz1 - vz2) * ((vz2.z() - m_zmax) / (vz1.z() - vz2.z()) + DBL_EPSILON);

	// ��������� ������ ����������� � ��������� ������ ����� ����������
	// � �� �������� ����������� ��������� ������ ����������� ���������� Z
	// (��� ��� ������� ��������������� �������� ����� Z),
	// �� ������� ���������� ����������� �� �������.

	// ���������! ��� ����������, ����� ������� ���������� � ����������� ���� ��� ������,
	// ��� ��������� ������������ ����� int.
	// ������� ��������� �� ������� �������� �� �� ��������, � �� ��������� � ������������.
	double rr = vz1.lengthXY() / (1.0 + (vz1.z() - m_iso) / m_sad);
	if (rr >= m_rmax)
	{
		// ���������� �� ����������� �������� ������
		double d = mcGeometry::getDistanceToConeOutside(geomVector3D(vz1.x(), vz1.y(), m_iso - vz1.z()), geomVector3D(v.x(), v.y(), -v.z()), m_rmax, m_sad);
		if (d == DBL_MAX)
			return;		// ���� ��, ���� ���� ������� �����������
		vz1 += v * (d + DBL_EPSILON);
		rr = m_rmax;
	}
	int ir = int(rr / m_rstep);
	if (ir >= m_nr) ir = m_nr - 1;

	// ��������� ������ �������
	double eddens = edep / (p1 - p0).length();
	int irnext = ir;

	while ((vz2 - vz1) * v > DBL_EPSILON)
	{
		ir = irnext;
		double r = ir * m_rstep;
		double d = DBL_MAX;

		// ���� ���� ����������� � ���������� ���������, �� �� �������� �� ���� � �� �� ����
		// � ��������� ����������� � ���������� � �� �������.
		if (ir > 0)
			d = mcGeometry::getDistanceToConeOutside(geomVector3D(vz1.x(), vz1.y(), m_iso - vz1.z()), geomVector3D(v.x(), v.y(), -v.z()), r, m_sad);

		// ���� ���������� ���������� �����, �� �������� �������� ������ ��������.
		if (d == DBL_MAX)
		{
			// ������������� �������� ������
			d = mcGeometry::getDistanceToConeInside(geomVector3D(vz1.x(), vz1.y(), m_iso - vz1.z()), geomVector3D(v.x(), v.y(), -v.z()), r + m_rstep, m_sad);
			irnext++;
		}
		else
			irnext--;

		if (d < 0)
			throw std::exception("mcScoreConicalRZ::ScoreLine: incorrect cones crossing");

		// ����� ����� � �������������� ������.
		// ������� ����� �� ���������� ������� ������, ��� ��������, ��� �� ���������� ������� ��������� ������ �������� ������.
		geomVector3D vz12 = d == DBL_MAX ? vz2 : vz1 + (v * (d + DBL_EPSILON));

		// ��������� �� ������ �� �� ������� ����� ����������� �����.
		if ((vz12 - vz2) * v > 0)
			vz12 = vz2;

		// ������� ���� � ������� ����������
		int iz = int((vz1.z() - m_zmin) / m_zstep);
		double z2 = m_zmin + m_zstep * (iz + 1);

		// ���� ����� ������� � ��� �� ����, �� � �� ���� ����������
		if (vz12.z() < z2)
			scoreEnergyInVoxel(iThread, ir, iz, eddens * (vz12 - vz1).length());
		else
		{
			// ���������� ���� �� Z
			while ((vz12 - vz1) * v > DBL_EPSILON)
			{
				geomVector3D pz2 = vz1 + (v * ((z2 - vz1.z()) / v.z()));
				if (pz2.z() > vz12.z())
					pz2 = vz12;

				scoreEnergyInVoxel(iThread, ir, iz, eddens * (pz2 - vz1).length());

				vz1 = pz2;
				iz++;
				z2 += m_zstep;
			}
		}
		vz1 = vz12;
	}
}

void mcScoreConicalRZ::scoreEnergyInVoxel(int iThread, int ir, int iz, double edep)
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreConicalRZ::scoreEnergyInVoxel: thread exceed container size");

	// ������� ��������� �������� ��� ����� ������� ��������, ����� ������ �� ������ ���������� �������������
	if (ir < 0)
	{
		std::strstream s;
		s << "mcScoreConicalRZ::scoreEnergyInVoxel: r index out of range" << std::endl;
		s << "iThread = \t" << iThread << std::endl;
		s << "ir = \t" << ir << std::endl;
		s << "iz = \t" << iz << std::endl;
		s << "edep = \t" << edep << std::endl;

		//throw std::exception(s.str());
		cout << s.str();
		return;
	}

	if (ir < m_nr && iz >= 0 && iz < m_nz) {
		m_M[iThread][ir*m_nz + iz] += edep;
		etotal_[iThread] += edep;
	}
}

void mcScoreConicalRZ::CE2D()
{
	if (dconverted_)
		throw std::exception("mcScoreConicalRZ::CE2D: Already converted to dose");
	dconverted_ = true;
	for (int iz = 0; iz < m_nz; iz++)
	{
		double z1 = m_zmin + iz * m_zstep - m_iso;
		double z2 = z1 + m_zstep;
		double fz = PI * m_zstep * (1 + (z2 + z1) / m_sad + (z2 * z2 + z1 * z2 + z1 * z1) / (m_sad * m_sad));
		for (int ir = 0; ir < m_nr; ir++)
		{
			double r1 = ir*m_rstep;
			double r2 = r1 + m_rstep;
			double vol = (r2*r2 - r1*r1) * fz;
			for (int i = 0; i < nThreads_; i++)
				m_M[i][ir*m_nz + iz] /= vol;
		}
	}
}

double mcScoreConicalRZ::Dose(int iThread, int ir, int iz) const //���� � ���������� ������
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreConicalRZ::Dose: thread exceed container size");
	else if (iz < 0 || iz >= m_nz || ir < 0 || ir >= m_nr)
		return 0;
	return m_M[iThread][ir*m_nz + iz];
}

double mcScoreConicalRZ::Dose(int ir, int iz) const //��������� ���� �� ������ (������� ���� ���������� �� ���-�� �������)
{
	double f = 0;
	if (iz >= 0 && iz < m_nz && ir >= 0 && ir < m_nr)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_M[i][ir*m_nz + iz];
	}
	return f;
}

double mcScoreConicalRZ::Sigma(int ir, int iz) const
{
	double f = Dose(ir, iz);// ��������� ���� �� ������ (������� ���� ����� �� ���-�� �������)
	if (f == 0)
		return 0;
	f /= nThreads_; // ������� ���� �� ������
	double sigma = 0;
	for (int register i = 0; i < nThreads_; i++)///���� ��� ������������ ����� �� ������
	{
		sigma += pow((m_M[i][ir*m_nz + iz] - f), 2);
	}
	sigma = sqrt(sigma / nThreads_);
	return sigma / f;
}

void mcScoreConicalRZ::dumpVRML(ostream& os) const
{
	os << "# mcScoreConicalRZ Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int ir, iz, it, count = 0;
	int da = 15; // ��� �� ���� 15 ��������
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

	for (iz = 0; iz <= m_nz; iz++)
	{
		double z = m_zmin + iz * m_zstep;
		// ��������������� �����
		for (ir = 1; ir <= m_nr; ir++) {
			double r = m_rstep * ir * (1.0 + (z - m_iso) / m_sad);
			for (it = 0; it < 360; it += da) {
				geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), z) * mttow;
				os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
				p = geomVector3D(r*sin(mPi*(it + da)), r*cos(mPi*(it + da)), z) * mttow;
				os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
				count++;
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

void mcScoreConicalRZ::dumpStatistic(ostream& os) const
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

	os << "NR\tNz\tRMAX\tZ1\tZ2\tZIso\tSad" << endl;
	os << nr << '\t' << nz << '\t' << rmax << '\t' << zmin << '\t' << zmax << '\t' << m_iso << '\t' << m_sad << endl << endl;

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

	// this code - for Sigma
	os << endl;
	os << "Sigma matrix" << endl;
	os << endl;
	for (ir = 0; ir < nr; ir++)
		os << '\t' << (double(ir) + 0.5) * rstep;
	os << endl;

	for (iz = 0; iz < nz; iz++)
	{
		os << zmin + (double(iz) + 0.5) * zstep;

		//"================Sigma==================";
		for (ir = 0; ir < nr; ir++)
			os << '\t' << Sigma(ir, iz);
		os << endl;
	}

	double sigmaMean = 0;
	int count = 0;
	for (iz = 4 * nz / 10; iz < 6 * nz / 10; iz++)
	{
		for (ir = 2 * nr / 10; ir < 4 * nr / 10; ir++)
		{
			sigmaMean += Sigma(ir, iz);
			count++;
		}
	}
	os << endl;
	os << "Mean sigma = \t" << sigmaMean / count << endl;
	os << endl;

	// ������ � ����������� �����.
	// ���� ���������� - �� ������� �� ���������� �� 1.5 �� �� ���� ���� 
	// �������������� ��������� � ������ ���� � ����������� ������.
	// �����, �� ������� ��������� SG ������� �������, ��� ��� ���� ������� � FanLine ��������� � 
	// � ������� �������������� �� ������� �� ������ ���� ������ ��������� �� ����������� ����� � ������.

	if (dconverted_)
	{
		// ��������� ������� ��������� �� threads, �� �������� ����� �������� ��������� �������.
		// ������� ��� � ����� ���������������.
		vector<vector<double>> nmatrix(nz, vector<double>(nr, 0));
		for (iz = 0; iz < nz; iz++)
			for (ir = 0; ir < nr; ir++)
				nmatrix[iz][ir] = Dose(ir, iz);

		auto smatrix = ProfileProcessor::SmoothFanRZ(nmatrix, nr, nz, rstep);

		os << endl;
		os << "Smooth dose matrix" << endl;
		os << endl;
		for (ir = 0; ir < nr; ir++)
			os << '\t' << (double(ir) + 0.5) * rstep;
		os << endl;

		for (iz = 0; iz < nz; iz++)
		{
			os << zmin + (double(iz) + 0.5) * zstep;
			for (ir = 0; ir < nr; ir++)
				os << '\t' << (*smatrix)[iz][ir];
			os << endl;
		}
	}
}
