#include "mcScore.h"

mcScore::mcScore(const char* module_name, int nThreads)
	:mcObj()
	, dconverted_(false)
	, density_(1.0)
	, transport_(nullptr)
{
	mcObj::setName("Default score");
	if (module_name != nullptr)
		strcpy_s(attached_module_, 64, module_name);
	nThreads_ = nThreads;
	etotal_ = new double[nThreads];
	for (int i = 0; i < nThreads; i++) etotal_[i] = 0;
}

mcScore::~mcScore()
{
	if (etotal_ != nullptr)
		delete[] etotal_;
}

double mcScore::etotal() const
{
	double etotal = 0;
	for (int i = 0; i < nThreads_; i++) etotal += etotal_[i];
	return etotal;
}

void mcScore::dumpVRML(ostream& os) const
{
	os << "# VRML dump for scoring is not implemented" << endl;
}

void mcScore::dumpStatistic(ostream& os) const
{
	os << *(const mcObj*)this << endl;
	os << "--------------------------------------------------------------------------" << endl;
	os << "Attached to module: \t" << getModuleName() << endl;
	os << "--------------------------------------------------------------------------" << endl;
	os << "Etotal = " << etotal() << endl;
	os << endl;
}
