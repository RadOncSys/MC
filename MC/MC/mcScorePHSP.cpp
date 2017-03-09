#include "mcScorePHSP.h"
#include "mcPhaseSpaceIO.h"

mcScorePHSP::mcScorePHSP(const char* module_name, const char* fname)
	:mcScore(module_name, 1)
{
	// Открыть файл для записи
	phsp_ = new mcPhaseSpaceIO();
	phsp_->open(fname, mcPhaseSpaceIO::WRITING);
}

mcScorePHSP::~mcScorePHSP(void)
{
	delete phsp_;
}

void mcScorePHSP::ScoreFluence(const mcParticle& particle)
{
	etotal_[0] += particle.ke * particle.weight;
	phsp_->write(particle);
}

void mcScorePHSP::dumpVRML(ostream& os) const
{
	os << "# VRML dump for mcScorePHSP scoring is not implemented" << endl;
}

void mcScorePHSP::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);
}
