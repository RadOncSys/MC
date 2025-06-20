#include "mcScoreFluence3D.h"
#include "mcGeometry.h"
#include "mctransport.h"

mcScoreFluence3D::mcScoreFluence3D(const char* module_name, int nThreads, 
	int nx, int ny, int nz, double psx, double psy, double psz, double z0)
	:mcScore(module_name, nThreads), nx_(nx), ny_(ny), nz_(nz), 
	psx_(psx), psy_(psy), psz_(psz), z0_(z0)
{
	int n = nx * ny * nz;
	M_ = new double* [nThreads];
	for (int i = 0; i < nThreads_; i++)
		M_[i] = new double[n];
}

mcScoreFluence3D::~mcScoreFluence3D()
{
	for (int i = 0; i < nThreads_; i++)
		delete [] M_[i];
	if (M_) delete[] M_;
}

void mcScoreFluence3D::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)
{
	return;
}

void mcScoreFluence3D::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);
}

void mcScoreFluence3D::dumpVRML(ostream& os)const
{
	os << "# Score: " << name_ << endl;
	if (transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

}
