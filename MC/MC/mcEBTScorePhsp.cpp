#include "mcEBTScorePhsp.h"
#include "mcThread.h"
#include <fstream>

mcEBTScorePhsp::mcEBTScorePhsp(const char* module_name, int nThreads, const char* outfname)
	:mcScore(module_name, nThreads), outfname_(outfname), npars_(7)
{
	phsp_.resize(nThreads);
}

mcEBTScorePhsp::~mcEBTScorePhsp()
{
}

void mcEBTScorePhsp::ScoreFluence(const mcParticle& particle)
{
	int iThread = particle.thread_->id();
	double edep = particle.ke * particle.weight;

	if (particle.t != MCP_PHOTON)
		return;

	phsp_[iThread].push_back(particle.ke);
	phsp_[iThread].push_back(particle.p.p_[0]);
	phsp_[iThread].push_back(particle.p.p_[1]);
	phsp_[iThread].push_back(particle.p.p_[2]);
	phsp_[iThread].push_back(particle.u.p_[0]);
	phsp_[iThread].push_back(particle.u.p_[1]);
	phsp_[iThread].push_back(particle.u.p_[2]);
}

void mcEBTScorePhsp::SavePHSP() const
{
	/*
	ofstream os(outfname_.c_str());
	if (os.fail())
		throw exception((std::string("mcEBTScorePhsp: Can't open file for writing particles - ") + outfname_).c_str());

	for (int it = 0; it < nThreads_; it++)
	{
		auto& pl = phsp_[it];
		for (int i = 0; i < pl.size();)
		{
			for (int j = 0; j < npars_; j++, i++)
			{
				os << pl[i];
				if (j < npars_ - 1)
					os << "\t";
				else
					os << std::endl;
			}
		}
	}
	*/

	int nParticleValues = 0;
	for (int it = 0; it < nThreads_; it++)
		nParticleValues += phsp_[it].size();

	int size = sizeof(int) + nParticleValues * sizeof(double);
	void* buffer = malloc(size);
	char* pbuffer = (char*)buffer;

	// Первое число - количество частиц в файле
	*(int*)pbuffer = int(nParticleValues / npars_);
	pbuffer += sizeof(int);
	
	for (int it = 0; it < nThreads_; it++)
	{
		int n = phsp_[it].size() * sizeof(double);
		memcpy(pbuffer, &phsp_[it][0], n);
		pbuffer += n;
	}

	FILE* file;
	if (fopen_s(&file, outfname_.c_str(), "wb") != 0)
		throw std::exception("mcEBTScorePhsp:: Cannot open model file for writing");
	fwrite(buffer, 1, size, file);
	fclose(file);

	free(buffer);
}

void mcEBTScorePhsp::dumpVRML(ostream& os) const
{
	os << "# mcEBTScorePhsp Score: " << name_ << endl;
	os << "# Not implemented!" << endl;

	SavePHSP();
}

void mcEBTScorePhsp::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);
}
