#include "mcScoreMatrixXY.h"
#include "mcGeometry.h"
#include "mctransport.h"
#include "ProfileProcessor.h"
//#include "../../Include/SimpleImage.h"

mcScoreMatrixXY::mcScoreMatrixXY(const char* module_name, int nThreads, int nx, int ny, double psx, double psy, double psz, double z0)
	:mcScore(module_name, nThreads), nx_(nx), ny_(ny), z0_(z0), psx_(psx), psy_(psy), psz_(psz)
{
	z1_ = z0_ + psz_;
	minx_ = -0.5 * nx_ * psx_;
	miny_ = -0.5 * ny_ * psy_;

	int len = nThreads * nx_ * ny_;
	MAll_ = new double[len];
	memset(MAll_, 0, len * sizeof(double));
	M_ = new double*[nThreads];
	M_[0] = MAll_;
	for (int i = 1; i < nThreads_; i++)
		M_[i] = M_[i - 1] + nx_ * ny_;
}

mcScoreMatrixXY::~mcScoreMatrixXY()
{
	if (MAll_) delete[] MAll_;
	if (M_) delete[] M_;
}

void mcScoreMatrixXY::ScorePoint(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p)
{
	if (p.z() < z0_ || p.z() >= z1_) return;
	if (p.x() <= minx_ || p.y() <= miny_)
		return;

	int ix = int((p.x() - minx_) / psx_);
	int iy = int((p.y() - miny_) / psy_);
	if (ix >= nx_ || iy >= ny_)
		return;

	M_[iThread][iy*nx_ + ix] += edep;
	etotal_[iThread] += edep;
}

void mcScoreMatrixXY::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)
{
	double x1 = p0.x(), y1 = p0.y(), z1 = p0.z();
	double x2 = p1.x(), y2 = p1.y(), z2 = p1.z();

	// Make program faster
	if ((z1 <= z0_ && z2 <= z0_) || (z1 >= z1_ && z2 >= z1_))
		return;
	if ((x1 <= minx_ && x2 <= minx_) || (y1 <= miny_ && y2 <= miny_))
		return;

	// Index limits
	int ix1 = (x1 < minx_) ? -1 : int((x1 - minx_) / psx_);
	int ix2 = (x2 < minx_) ? -1 : int((x2 - minx_) / psx_);
	int iy1 = (z1 < miny_) ? -1 : int((y1 - miny_) / psy_);
	int iy2 = (z2 < miny_) ? -1 : int((y2 - miny_) / psy_);
	if ((ix1 >= nx_ && ix2 >= nx_) || (iy1 >= ny_ && iy2 >= ny_))
		return;

	// When path is completely within one voxel
	if (ix1 == ix2 && iy1 == iy2 &&
		z1 >= z0_ && z2 >= z0_ && z1 <= z1_ && z2 <= z1_)
	{
		M_[iThread][iy1*nx_ + ix1] += edep;
		etotal_[iThread] += edep;
		return;
	}

	ix1 = (ix1 < 0) ? 0 : (ix1 >= nx_) ? nx_ - 1 : ix1;
	ix2 = (ix2 < 0) ? 0 : (ix2 >= nx_) ? nx_ - 1 : ix2;
	iy1 = (iy1 < 0) ? 0 : (iy1 >= ny_) ? ny_ - 1 : iy1;
	iy2 = (iy2 < 0) ? 0 : (iy2 >= ny_) ? ny_ - 1 : iy2;

	int i;
	if (ix1 > ix2) { i = ix1; ix1 = ix2; ix2 = i; }
	if (iy1 > iy2) { i = iy1; iy1 = iy2; iy2 = i; }

	// Energy lose per unit path.
	double f = edep / (p1 - p0).length();

	// Go trough voxels
	double zv1 = z0_, zv2 = z1_;
	for (int iy = iy1; iy <= iy2; iy++)
	{
		double yv1 = miny_ + psy_ * iy, yv2 = yv1 + psy_;
		for (int ix = ix1; ix <= ix2; ix++)
		{
			double xv1 = minx_ + psx_ * ix, xv2 = xv1 + psx_;
			double trl = mcGeometry::TrLenInVoxel(x1, y1, z1, x2, y2, z2, xv1, yv1, zv1, xv2, yv2, zv2);
			M_[iThread][iy*nx_ + ix] += f * trl;
			etotal_[iThread] += f * trl;
		}
	}
	return;
}

double mcScoreMatrixXY::Dose(int iThread, int ix, int iy) const
{
	if (iThread >= nThreads_)
		throw std::exception("mcScoreMatrixXY::Dose: thread exceed container size");
	else if (iy < 0 || iy >= ny_ || ix < 0 || ix >= nx_)
		return 0;
	return M_[iThread][iy*nx_ + ix];
}

double mcScoreMatrixXY::Dose(int ix, int iy) const
{
	double f = 0;
	if (iy >= 0 && iy < ny_ && ix >= 0 && ix < nx_)
	{
		for (int i = 0; i < nThreads_; i++)
			f += M_[i][iy*nx_ + ix];
	}
	return f;
}

void mcScoreMatrixXY::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	int i, j;
	if (ny_ == 1)
	{
		double width, penumbra, dindex;
		std::vector<double> x, d, dd;
		for (i = 0; i < nx_; i++)
		{
			x.push_back(minx_ + (double(i) + 0.5) * psx_);
			d.push_back(Dose(i, 0));
		}
		ProfileProcessor::ProfileParameters(x, d, dd, width, penumbra, dindex);

		os << endl;
		os << "Profile parameters" << endl;
		os << "------------------" << endl;
		os << "Width = \t" << width << endl;
		os << "Penumbra = \t" << penumbra << endl;
		os << "DIndex = \t" << dindex << endl;

		os << nx_ << '\t' << ny_ << '\t' << psx_ << '\t' << psy_ << z0_ << endl << endl;
		for (i = 0; i < nx_; i++)
			os << x[i] << '\t' << d[i] << '\t' << dd[i] << endl;
	}
	else
	{
		os << "Energy distribution in XY plane:" << endl;

		os << "NX\tNY\tPSX\tPSY\tZ0" << endl;
		os << nx_ << '\t' << ny_ << '\t' << psx_ << '\t' << psy_ << z0_ << endl << endl;

		for (i = 0; i < nx_; i++)
			os << '\t' << minx_ + (double(i) + 0.5) * psx_;
		os << endl;

		for (j = 0; j < ny_; j++) {
			os << miny_ + (double(j) + 0.5) * psy_;
			for (i = 0; i < nx_; i++)
				os << '\t' << Dose(i, j);
			os << endl;
		}
	}
}

void mcScoreMatrixXY::dumpVRML(ostream& os)const
{
	os << "# Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	unsigned i, j, k, ni = nx_ + 1, nj = ny_ + 1, nk = 2;

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

	// Грани XY (четные точки на Zmin, нечетные - Zmax)
	unsigned ncount = 0;
	for (j = 0; j < nj; j++) {
		if (j != 0 && j != nj - 1) continue;	// Пропускаем внутренние линии. Рисуем только периметр
		double y = miny_ + j * psy_;
		for (i = 0; i < ni; i++) {
			if (i != 0 && i != ni - 1) continue;
			double x = minx_ + i * psx_;
			geomVector3D p = geomVector3D(x, y, z0_) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(x, y, z1_) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			ncount++;
		}
	}
	// Грани XZ
	for (k = 0; k < nk; k++) {
		double z = z0_ + k * psz_;
		for (i = 0; i < ni; i++) {
			if (i != 0 && i != ni - 1) continue;
			double x = minx_ + i * psx_;
			geomVector3D p = geomVector3D(x, miny_, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(x, miny_ + psy_*ny_, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			ncount++;
		}
	}
	// Грани YZ
	for (k = 0; k < nk; k++) {
		double z = z0_ + k * psz_;
		for (j = 0; j < nj; j++) {
			if (j != 0 && j != nj - 1) continue;
			double y = miny_ + j * psy_;
			geomVector3D p = geomVector3D(minx_, y, z) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(minx_ + psx_*nx_, y, z) * mttow;
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
