#include "mcGaussSource.h"
#include "mcDefs.h"
#include "mcSamplers.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

mcGaussSource::mcGaussSource(const char* name, int nThreads, mc_particle_t type, double ke, double spreadke, double z, double sigmax, double sigmay, double sigmathetax, double sigmathetay)
	:mcSource(name, nThreads)
	, type_(type)
	, ke_(ke)
	, z_(z)
	, spreadke_(spreadke) // ğàçáğîñ êèíåòè÷åñêîé ıíåğãèè â ïğîöåíòàõ
	, sigmax_(sigmax* sqrt(0.5))
	, sigmay_(sigmay* sqrt(0.5))
	, sigmathetax_(sigmathetax* sqrt(0.5))
	,sigmathetay_(sigmathetay* sqrt(0.5))

{
	q_ = (type_ == mc_particle_t::MCP_NEGATRON) ? -1 : (type_ == mc_particle_t::MCP_POSITRON || type_ == mc_particle_t::MCP_PROTON) ? 1 : 0;
}

mcGaussSource::~mcGaussSource(void)
{
}

void mcGaussSource::sample(mcParticle& p, mcThread* thread)
{
	
	const double s = sqrt(0.5);
	mcRng& rng = thread->rng();

	p.t = type_;
	p.q = q_;
	p.ke = ke_+ke_*spreadke_*0.01*(-1 + 2  * rng.rnd());

	if (sigmax_ == 0 && sigmay_ == 0)
		p.p.set(0, 0, z_);
	else
	{
		double x = s * sigmax_ * mcSamplers::SampleGauss2D(rng.rnd()); // Ğîçûãğûø ìåñòà ğîæäåíèÿ ÷àñòèö â ıëëèïìå íà îñè Õ
		double y = s * sigmay_ * mcSamplers::SampleGauss2D(rng.rnd()); // Ğîçûãğûø  ìåñòà ğîæäåíèÿ ÷àñòèö â ıëëèïñå ïî îñè y
		p.p.set(x, y, z_);
		//p.p.set(0, 0, z_);
	}
	p.plast = p.p;

	
		double phix = s * sigmathetax_ * mcSamplers::SampleGauss2D(rng.rnd());
		double phiy = s * sigmathetay_ * mcSamplers::SampleGauss2D(rng.rnd());
		p.u.set(sin(phix), sin(phiy) , sqrt(1- sin(phix)*sin(phix)- sin(phiy)*sin(phiy)));
		//p.u.set(0, 0, 1);
	

	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += ke_;
}

