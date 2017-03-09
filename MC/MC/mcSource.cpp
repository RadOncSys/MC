#include "mcsource.h"
#include "mcScoreTrack.h"

mcSource::mcSource()
	:mcObj()
	, trackScore_(nullptr), nThreads_(0), etotal_(nullptr), isGamma_(false)
{
	red_ = 1.0; green_ = 0; blue_ = 0; transparancy_ = 0.2;
}

mcSource::mcSource(const char* name, int nThreads)
	:mcObj()
	, trackScore_(nullptr)
	, isGamma_(false)
{
	red_ = 1.0; green_ = 0; blue_ = 0; transparancy_ = 0.2;
	setName(name);
	nThreads_ = nThreads;
	etotal_ = new double[nThreads];
	for (int i = 0; i < nThreads; i++) etotal_[i] = 0;
}

mcSource::~mcSource()
{
	if (etotal_ != nullptr) delete[] etotal_;
	if (trackScore_) delete trackScore_;
}

void mcSource::setScoreTrack(double R, double Z1, double Z2, double EMIN, bool doPhotons, bool doElectrons, bool doPositrons)
{
	trackScore_ = new mcScoreTrack(nThreads_, R, Z1, Z2, EMIN, doPhotons, doElectrons, doPositrons);
}

void mcSource::setModuleName(const char* s)
{
	strcpy_s(attached_module_, 64, s);
}

double mcSource::etotal() const
{
	double etotal = 0;
	for (int i = 0; i < nThreads_; i++) etotal += etotal_[i];
	return etotal;
}

void mcSource::dumpVRML(ostream& os) const
{
	if (trackScore_)
		trackScore_->dumpVRML(os);
}
