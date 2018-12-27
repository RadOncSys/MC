#include "mcParticle.h"

mcParticle::mcParticle(void)
	:t(MCP_PHOTON)
	, q(0)
	, ke(0)
	, dnear(0)
	, weight(1)
	, mfps(0)
	, regDensityRatio(1)
	, transport_(nullptr)
	, transportNearest_(nullptr)
	, trackScore_(nullptr)
	, thread_(nullptr)
	, regionBirth(0)
	, regionFlags(0)
	, exitSurface_(temb_shit_t::Undefined)
{
}

mcParticle::mcParticle(mc_particle_t pt, int pq, double pke, const geomVector3D& pp, const geomVector3D& pu)
	:t(pt)
	, q(pq)
	, ke(pke)
	, p(pp)
	, u(pu)
	, dnear(0)
	, weight(1)
	, mfps(0)
	, regDensityRatio(1)
	, transport_(nullptr)
	, transportNearest_(nullptr)
	, trackScore_(nullptr)
	, thread_(nullptr)
	, regionBirth(0)
	, regionFlags(0)
	, exitSurface_(temb_shit_t::Undefined)
{
}

mcParticle::mcParticle(const mcParticle& p)
	:t(p.t)
	, q(p.q)
	, ke(p.ke)
	, p(p.p)
	, plast(p.plast)
	, u(p.u)
	, dnear(p.dnear)
	, region(p.region)
	, weight(p.weight)
	, mfps(p.mfps)
	, regDensityRatio(p.regDensityRatio)
	, transport_(p.transport_)
	, transportNearest_(p.transportNearest_)
	, trackScore_(p.trackScore_)
	, thread_(p.thread_)
	, regionBirth(p.regionBirth)
	, regionFlags(p.regionFlags)
	, exitSurface_(p.exitSurface_)
{
}

mcParticle::~mcParticle(void)
{
}
