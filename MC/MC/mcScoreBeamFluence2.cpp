#include "mcScoreBeamFluence2.h"
#include "mcParticle.h"
#include "mcThread.h"
#include "mcTransport.h"

mcScoreBeamFluence2::mcScoreBeamFluence2(const char* module_name, int nThreads, int nr, double rmax, int nr_s, double rmax_s, double H, int nz, double d_iz)
	:mcScore(module_name, nThreads), h_(H), nz_(nz), d_iz_(d_iz)
{
	m_nr = nr;
	m_rmax = rmax;
	m_rstep = m_rmax / m_nr;

	int len = nThreads * m_nr;
	m_intencity_all = new double[len];
	m_intencity_count_all = new double[len];
	memset(m_intencity_all, 0, len * sizeof(double));
	memset(m_intencity_count_all, 0, len * sizeof(double));
	m_intencity = new double*[nThreads];
	m_intencity_count = new double*[nThreads];
	m_intencity[0] = m_intencity_all;
	m_intencity_count[0] = m_intencity_count_all;
	for (int i = 1; i < nThreads_; i++)
	{
		m_intencity[i] = m_intencity[i - 1] + m_nr;
		m_intencity_count[i] = m_intencity_count[i - 1] + m_nr;
	}

	m_nr_s = nr_s;
	m_rmax_s = rmax_s;
	m_rstep_s = m_rmax_s / m_nr_s;

	// Трехмерный массив скоринга фотографии источника на серии уровней
	m_intencity_s = new double**[nThreads];
	m_intencity_s[0] = new double*[nThreads * nz_];
	int len_s = nThreads * nz_ * m_nr_s;
	m_intencity_s[0][0] = new double[len_s];
	if (!m_intencity_s || !m_intencity_s[0] || !m_intencity_s[0][0])
		throw std::exception("mcScoreBeamFluence2: memory error");
	memset(m_intencity_s[0][0], 0, len_s * sizeof(double));

	int it, iz;
	for (it = 1; it < nThreads; it++)
		m_intencity_s[it] = m_intencity_s[it - 1] + nz_;
	for (iz = 1; iz < nz; iz++)
		m_intencity_s[0][iz] = m_intencity_s[0][iz - 1] + m_nr_s;
	for (it = 1; it < nThreads; it++)
	{
		m_intencity_s[it][0] = m_intencity_s[it - 1][0] + nz * m_nr_s;
		for (iz = 1; iz < nz; iz++) m_intencity_s[it][iz] = m_intencity_s[it][iz - 1] + m_nr_s;
	}
}

mcScoreBeamFluence2::~mcScoreBeamFluence2()
{
	if (m_intencity_all) delete[] m_intencity_all;
	if (m_intencity_count_all) delete[] m_intencity_count_all;
	if (m_intencity) delete[] m_intencity;
	if (m_intencity_count) delete[] m_intencity_count;

	if (m_intencity_s) {
		delete[] m_intencity_s[0][0];
		delete[] m_intencity_s[0];
		delete[] m_intencity_s;
	}
}

void mcScoreBeamFluence2::ScoreFluence(const mcParticle& particle)
{
	int iThread = particle.thread_->id();
	double edep = particle.ke * particle.weight;
	etotal_[iThread] += edep;

	// HACK !!!	GG 20140703
	if (particle.t != MCP_PHOTON)
		return;

	// !!! Scoring вызывается до перемещения частицы на поверхность объекта к которому привязана частица.
	// Нужно ее перенести здесь на поверхность.
	geomVector3D p = particle.p + (particle.u * (-particle.p.z() / particle.u.z()));

	double r = p.lengthXY();
	if (r >= m_rmax)
		return;
	int ir = int(p.lengthXY() / m_rstep);
	if (ir < m_nr)
	{
		m_intencity[iThread][ir] += edep;
		m_intencity_count[iThread][ir] += particle.weight;
	}

	// Проектируем частицу на плоскость источника
	for (int iz = 0; iz < nz_; iz++)
	{
		geomVector3D pc = p + (particle.u * ((-h_ + iz*d_iz_) / particle.u.z()));

		r = pc.lengthXY();
		if (r >= m_rmax_s)
			return;
		ir = int(pc.lengthXY() / m_rstep_s);
		if (ir < m_nr_s)
			m_intencity_s[iThread][iz][ir] += edep;
	}
}

double mcScoreBeamFluence2::Intencity(int iThread, int ir) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluence2::Intencity: thread exceed container size");
	else if (ir < 0 || ir >= m_nr)
		return 0;
	return m_intencity[iThread][ir] / (PI * m_rstep * m_rstep * (2 * ir + 1));
}

double mcScoreBeamFluence2::MeanEnergy(int iThread, int ir) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluence2::Intencity: thread exceed container size");
	else if (ir < 0 || ir >= m_nr)
		return 0;
	return (m_intencity_count[iThread][ir] > 0) ? m_intencity[iThread][ir] / m_intencity_count[iThread][ir] : 0;
}

double mcScoreBeamFluence2::Intencity_s(int iThread, int iz, int ir) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreBeamFluence2::Intencity_s: thread exceed container size");
	else if (ir < 0 || ir >= m_nr_s)
		return 0;
	else
	{
		return m_intencity_s[iThread][iz][ir] / (PI * m_rstep_s * m_rstep_s * (2 * ir + 1));
	}
}

double mcScoreBeamFluence2::Intencity(int ir) const
{
	double f = 0;
	if (ir >= 0 && ir < m_nr)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_intencity[i][ir];
		f /= PI * m_rstep * m_rstep * (2 * ir + 1);
	}
	return f;
}

double mcScoreBeamFluence2::MeanEnergy(int ir) const
{
	double f = 0, fcount = 0;
	if (ir >= 0 && ir < m_nr)
	{
		for (int i = 0; i < nThreads_; i++)
		{
			f += m_intencity[i][ir];
			fcount += m_intencity_count[i][ir];
		}
	}
	return fcount > 0 ? f / fcount : 0;
}

double mcScoreBeamFluence2::Intencity_s(int iz, int ir) const
{
	double f = 0;
	if (ir >= 0 && ir < m_nr_s && iz >= 0 && iz < nz_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_intencity_s[i][iz][ir];
		f /= PI * m_rstep_s * m_rstep_s * (2 * ir + 1);
	}
	return f;
}

double mcScoreBeamFluence2::Value_s(int iz, int ir) const
{
	double f = 0;
	if (ir >= 0 && ir < m_nr_s && iz >= 0 && iz < nz_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += m_intencity_s[i][iz][ir];
	}
	return f;
}

void mcScoreBeamFluence2::dumpVRML(ostream& os) const
{
	os << "# mcScoreBeamFluence2 Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int ir, it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;

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
	for (ir = 1; ir <= m_nr; ir++) {
		double r = m_rstep * ir;
		for (it = 0; it < 360; it += da) {
			geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), 0) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r*sin(mPi*(it + da)), r*cos(mPi*(it + da)), 0) * mttow;
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

void mcScoreBeamFluence2::dumpStatistic(ostream& os) const
{
	int ir, iz;

	mcScore::dumpStatistic(os);

	os << endl;
	os << "Mean energy profile" << endl;
	os << endl;

	for (ir = 0; ir < m_nr; ir++)
		os << '\t' << (double(ir) + 0.5) * m_rstep;
	os << endl;

	for (ir = 0; ir < m_nr; ir++)
		os << '\t' << MeanEnergy(ir);
	os << endl;

	os << endl;
	os << "Fluence intencity profile" << endl;
	os << endl;

	for (ir = 0; ir < m_nr; ir++)
		os << '\t' << (double(ir) + 0.5) * m_rstep;
	os << endl;

	for (ir = 0; ir < m_nr; ir++)
		os << '\t' << Intencity(ir);
	os << endl;

	os << endl;
	os << "Source intencity profile" << endl;
	os << endl;

	for (ir = 0; ir < m_nr_s; ir++)
		os << '\t' << (double(ir) + 0.5) * m_rstep_s;
	os << endl;

	for (iz = 0; iz < nz_; iz++)
	{
		os << (-h_ + iz*d_iz_);
		for (ir = 0; ir < m_nr_s; ir++)
			os << '\t' << Intencity_s(iz, ir);
		os << endl;
	}

	os << endl;
	os << "Comulative source histogram" << endl;
	os << endl;

	for (ir = 0; ir <= m_nr_s; ir++)
		os << '\t' << double(ir) * m_rstep_s;
	os << endl;

	for (iz = 0; iz < nz_; iz++)
	{
		double sum = 0;
		for (ir = 0; ir < m_nr_s; ir++)
			sum += Value_s(iz, ir);

		double count = 0;
		os << (-h_ + iz*d_iz_);
		os << '\t' << count;
		for (ir = 0; ir < m_nr_s; ir++)
		{
			count += Value_s(iz, ir);
			os << '\t' << count / sum;
		}
		os << endl;
	}
}
