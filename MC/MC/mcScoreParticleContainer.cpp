#include "mcScoreParticleContainer.h"
#include "mcThread.h"

mcScoreParticleContainer::mcScoreParticleContainer(const char* module_name, int nThreads)
	:mcScore(module_name, nThreads)
	, ptypeFilter_(MCP_NTYPES)
{
	particles_.resize(nThreads);
	etotal_other_ = new double[nThreads];
	for (int i = 0; i < nThreads; i++) etotal_other_[i] = 0;
}

mcScoreParticleContainer::~mcScoreParticleContainer(void)
{
	if (etotal_other_ != nullptr) delete[] etotal_other_;
}

void mcScoreParticleContainer::ScoreFluence(const mcParticle& particle)
{
	int iThread = particle.thread_->id();
	double edep = particle.ke * particle.weight;

	if (particle.t != ptypeFilter_)
	{
		etotal_other_[iThread] += edep;
		return;
	}
	else
		etotal_[iThread] += edep;

	PlaneParticleRecord pr;
	pr.x = (float)particle.p.x();
	pr.y = (float)particle.p.y();
	pr.z = (float)particle.p.z();
	pr.ux = (float)particle.u.x();
	pr.uy = (float)particle.u.y();
	pr.uz = (float)particle.u.z();
	pr.e = (float)particle.ke;
	pr.r = (float)particle.p.length();

	double a = particle.u.z();
	a = a < 0 ? 0 : a > 1 ? 1 : a;
	pr.a = (float)(acos(a) * 180 / PI);

	particles_[iThread].push_back(pr);
}

double mcScoreParticleContainer::etotal_other() const
{
	double etotal = 0;
	for (int i = 0; i < nThreads_; i++) etotal += etotal_other_[i];
	return etotal;
}

void mcScoreParticleContainer::dumpVRML(ostream& os) const
{
	os << "# mcScoreParticleContainer Score: " << name_ << endl;
	os << "# Not implemented!" << endl;
}

void mcScoreParticleContainer::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);

	os << "Etotal Others = " << etotal_other() << endl;

	// Расчетные параметры: распределение средней энергии по углу и распределение среднего угла по энергии.
	// Разбиваем частицы с шагом 1 градус и 0.01 МэВ
	float de = 0.1f;
	int eidxmax = 0;
	vector<float> emeans(90, 0.f);
	vector<float> ameans(size_t(25.f / de), 0.f);
	vector<size_t> ecounts(emeans.size(), 0);
	vector<size_t> acounts(ameans.size(), 0);

	int it;
	size_t i;
	for (it = 0; it < nThreads_; it++)
	{
		const vector<PlaneParticleRecord>& pl = particles_[it];
		for (i = 0; i < pl.size(); i++)
		{
			const PlaneParticleRecord& p = pl[i];
			size_t eidx = size_t(p.e / de);
			size_t aidx = size_t(p.a);

			if (aidx < emeans.size())
			{
				emeans[aidx] += p.e;
				ecounts[aidx]++;
			}

			if (eidx < ameans.size())
			{
				if ((int)eidx > eidxmax) eidxmax = (int)eidx;
				ameans[eidx] += p.a;
				acounts[eidx]++;
			}
		}
	}

	// Вывод средних значений

	os << "Mean energy distribution" << endl;
	os << "------------------------" << endl;
	os << "Angle\tEnergy" << endl;
	for (i = 0; i < emeans.size(); i++)
		os << (i + 0.5f) << "\t" << (ecounts[i] > 0 ? emeans[i] / ecounts[i] : 0) << endl;
	os << endl;

	os << "Mean angle distribution" << endl;
	os << "-----------------------" << endl;
	os << "Energy\tAngle" << endl;
	for (i = 0; i <= (size_t)eidxmax; i++)
		os << (de * (i + 0.5f)) << "\t" << (acounts[i] > 0 ? ameans[i] / acounts[i] : 0) << endl;
	os << endl;

	// Вывод частиц

	os << "List of" << (ptypeFilter_ == mc_particle_t::MCP_NEGATRON ? "electron" :
		ptypeFilter_ == mc_particle_t::MCP_POSITRON ? "positron" : 
		ptypeFilter_ == mc_particle_t::MCP_PROTON ? "proton" : 
		ptypeFilter_ == mc_particle_t::MCP_NEUTRON ? "neutron" : "photon")
		<< "s on the score surface:" << endl;
	os << endl;
	os << "X\tY\tZ\tUx\tUy\tUz\tE\tR\tAngle (degrees)" << endl;

	for (it = 0; it < nThreads_; it++)
	{
		const vector<PlaneParticleRecord>& pl = particles_[it];
		for (i = 0; i < pl.size(); i++)
		{
			const PlaneParticleRecord& p = pl[i];
			os << p.x << "\t" << p.y << "\t" << p.z << "\t"
				<< p.ux << "\t" << p.uy << "\t" << p.uz << "\t"
				<< p.e << "\t" << p.r << "\t" << p.a << endl;
		}
	}
	os << endl;
}
