#include "mcPhysicsNeutron.h"
#include "mcMediumNeutron.h"
#include "mcPhysicsCommon.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcDefs.h"
#include <float.h>

mcPhysicsNeutron::mcPhysicsNeutron(void)
{
}

mcPhysicsNeutron::~mcPhysicsNeutron(void)
{
}

bool mcPhysicsNeutron::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	if (p->ke <= ((const mcMediumNeutron&)med).transCutoff_neutron)
	{
		edep = p->ke;
		DiscardParticle(p);
		return true;
	}
	else
		return false;
}

double mcPhysicsNeutron::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	const mcMediumNeutron& m = (const mcMediumNeutron&)med;
	double logKE = ke;
	int iLogKE = int(ke);
	double sigma = (m.sigma0_proto[iLogKE] + logKE * m.sigma1_proto[iLogKE]) * dens;
	return (sigma > 0.0) ? 1 / sigma : DBL_MAX;
}

double mcPhysicsNeutron::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	static const double minBetaSquared = 1.0e-8;
	static const double twoOverEsubS = 0.094315;

	const mcMediumNeutron& m = (const mcMediumNeutron&)med;
	double e_dep = 0;

	double logKE = p->ke;
	int iLogKE = int(logKE);
	double dedx = p->regDensityRatio * (m.dedx0_proto[iLogKE] + logKE * m.dedx1_proto[iLogKE]);

	e_dep = 0.01 * p->ke;
	double mE = p->ke - e_dep / 2.0;
	logKE = mE;
	iLogKE = int(logKE);
	dedx = p->regDensityRatio * (m.dedx0_proto[iLogKE] + logKE * m.dedx1_proto[iLogKE]);
	step = MIN(step, e_dep / dedx);
	e_dep = step * dedx;

	double rel = (p->ke + PMASS) * PMASS_1;
	double sig2 = p->regDensityRatio * m.dEdxStragglingGaussVarianceConstPart_ * step * (1 + SQUARE(rel));

	double r1, r2;
	GaussStandardRnd_by_Marsaglia(p->thread_->rng(), r1, r2);
	e_dep = e_dep - sqrt(sig2) * r1;
	e_dep = (e_dep < 0) ? 0 : (e_dep > p->ke) ? p->ke : e_dep;

	p->p += p->u * step;

	double tE = p->ke + PMASS;

	double l = p->regDensityRatio * step / m.radLength;
	double a1_2 = (1 + 0.038 * log(l)) / (betasq(PMASS, tE) * tE);
	double th0 = sqrt(184.96 * l * SQUARE(a1_2));
	double theta = fabs(r2 * th0);
	double phi = p->thread_->rng().rnd() * TWOPI;
	ChangeDirection(cos(theta), sin(theta), cos(phi), sin(phi), p->u);

	// p->mfps = 0;
	p->ke -= e_dep;
	return e_dep;
}

double mcPhysicsNeutron::DoInterruction(mcParticle* p, const mcMedium* med) const
{
	p->ke /= 3.0;
	return 2 * p->ke * p->weight;
}
