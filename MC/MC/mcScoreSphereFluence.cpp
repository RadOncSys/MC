#include "mcScoreSphereFluence.h"
#include "mcThread.h"

mcScoreSphereFluence::mcScoreSphereFluence(const char* module_name, int nThreads)
	:mcScore(module_name, nThreads)
{
	particles_.resize(nThreads);
	etotal_other_ = new double[nThreads];
	for (int i = 0; i < nThreads; i++) etotal_other_[i] = 0;
}

mcScoreSphereFluence::~mcScoreSphereFluence()
{
	if (etotal_other_ != nullptr) delete[] etotal_other_;
}

void mcScoreSphereFluence::ScoreFluence(const mcParticle& particle)
{
	int iThread = particle.thread_->id();
	double edep = particle.ke * particle.weight;

	if (particle.t != MCP_PHOTON)
	{
		etotal_other_[iThread] += edep;
		return;
	}
	else
		etotal_[iThread] += edep;

	SFParticleRecord pr;
	pr.x = (float)particle.p.x();
	pr.y = (float)particle.p.y();
	pr.z = (float)particle.p.z();
	pr.ux = (float)particle.u.x();
	pr.uy = (float)particle.u.y();
	pr.uz = (float)particle.u.z();
	pr.e = (float)particle.ke;
	pr.r = (float)particle.p.length();

	geomVector3D punit(particle.p);
	punit.normalize();
	double a = punit * particle.u;
	a = a < 0 ? 0 : a > 1 ? 1 : a;
	pr.a = (float)(acos(a) * 180 / PI);

	particles_[iThread].push_back(pr);
}

double mcScoreSphereFluence::etotal_other() const
{
	double etotal = 0;
	for (int i = 0; i < nThreads_; i++) etotal += etotal_other_[i];
	return etotal;
}

void mcScoreSphereFluence::dumpVRML(ostream& os) const
{
	os << "# mcScoreSphereFluence Score: " << name_ << endl;
	os << "# Not implemented!" << endl;
}

void mcScoreSphereFluence::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	os << "Etotal Others = " << etotal_other() << endl;

	os << "List of photons on the score surface:" << endl;
	os << endl;
	os << "X\tY\tZ\tUx\tUy\tUz\tE\tR\tAngle (degrees)" << endl;

	for (int it = 0; it < nThreads_; it++)
	{
		const vector<SFParticleRecord>& pl = particles_[it];
		for (size_t i = 0; i < pl.size(); i++)
		{
			const SFParticleRecord& p = pl[i];
			os << p.x << "\t" << p.y << "\t" << p.z << "\t"
				<< p.ux << "\t" << p.uy << "\t" << p.uz << "\t"
				<< p.e << "\t" << p.r << "\t" << p.a << endl;
		}
	}
	os << endl;
}
