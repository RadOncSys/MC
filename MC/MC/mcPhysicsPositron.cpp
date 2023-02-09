#include "mcPhysicsPositron.h"
#include "mcMediumXE.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcTransport.h"
#include "mcDefs.h"
#include <float.h>

mcPhysicsPositron::mcPhysicsPositron(void)
{
}

mcPhysicsPositron::~mcPhysicsPositron(void)
{
}

bool mcPhysicsPositron::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	if (p->ke <= ((const mcMediumXE&)med).transCutoff_elec || p->ke <= p->transport_->transCutoff_elec)
	{
		edep = p->ke;
		AnnihilateAtRest(p);

		// Частицу не уничтожаем, хотя говорим, что да.
		// Да означает дальнейший учет оставшейся кинетической энергии позитрона.
		// Не уничтожаем, потому что ее место занимает один из фотонов в функции аннигиляции,
		// что эквивалентно уничтожению и созданию в стэке двух новых частиц.
		return true;
	}
	else
		return false;
}

double mcPhysicsPositron::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	const mcMediumXE& m = (const mcMediumXE&)med;
	double logKE = log(ke);
	int iLogKE = int(m.iLogKE0_elec + logKE * m.iLogKE1_elec);
	double sigma = (m.sigma0_posi[iLogKE] + logKE * m.sigma1_posi[iLogKE])* dens;
	return (sigma > 0.0) ? 1 / sigma : DBL_MAX;
}

double mcPhysicsPositron::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	static const double minBetaSquared = 1.0e-8;
	static const double twoOverEsubS = 0.094315;

	const mcMediumXE& m = (const mcMediumXE&)med;
	double e_dep;

	double logKE = log(p->ke);
	int iLogKE = (int)(m.iLogKE0_elec + logKE * m.iLogKE1_elec);
	double dedx = p->regDensityRatio * (m.dedx0_posi[iLogKE] + logKE * m.dedx1_posi[iLogKE]);

	// Ситуация, когда частица заведомо не выйдет из области.
	double range = p->ke / dedx;
	if (range < step && range < p->dnear)
	{
		step = range;
		e_dep = p->ke;
		p->p += p->u * step;
		p->ke = 0;
		p->dnear -= step;
		return e_dep;
	}

	// Compute some useful kinematic parameters
	double eTotal = p->ke + EMASS;
	double betaSquared = p->ke * (p->ke + TWICE_EMASS) / SQUARE(eTotal);
	betaSquared = MAX(betaSquared, minBetaSquared);
	double tscat = twoOverEsubS * betaSquared * eTotal;
	tscat = m.radLength * SQUARE(tscat) / p->regDensityRatio;

	// Реальный шаг выбирается как минимальный из шага до дискретного события, 
	// шага непрерывных потерь и расстояние, на котором влияние рассеяния не станвится критичным для модели.
	double stepSize = (m.stepSize0_posi[iLogKE] + logKE * m.stepSize1_posi[iLogKE]) / p->regDensityRatio;
	step = MIN(stepSize, MIN(step, 0.3 * tscat));

	// Длина трека больше расстояния между точками и связана рассеивающей поспособностью.
	double pathLength = ReducedPathLength(step, tscat);

	// Потери в треке расчитываются в два шага, сначала для dedx текущей энергии, затем средней.
	e_dep = dedx * pathLength;
	e_dep = MIN(e_dep, p->ke);

	logKE = log(p->ke - 0.5 * e_dep);
	iLogKE = int(m.iLogKE0_elec + logKE * m.iLogKE1_elec);
	if (iLogKE < 0) iLogKE = 0;  // вблизи энергии поглощения возможна проблема индексов
	dedx = p->regDensityRatio * (m.dedx0_posi[iLogKE] + logKE * m.dedx1_posi[iLogKE]);
	e_dep = dedx * pathLength;
	e_dep = MIN(e_dep, p->ke);
	p->p += p->u * step;

	// Change direction in accordance with the Moliere theory
	MoliereScatter(p->thread_->rng(),
		m.scaledLittleB, m.chi_cc, p->u,
		p->ke + EMASS, betaSquared,
		pathLength * p->regDensityRatio);

	p->plast = p->p;
	// p->mfps = 0;
	p->ke -= e_dep;
	p->dnear -= step;
	return e_dep;
}

double mcPhysicsPositron::DoInterruction(mcParticle* p, const mcMedium* med) const
{
	const mcMediumXE* m = (const mcMediumXE*)med;
	double e_dep = 0;

	double logKE = log(p->ke);
	int iLogKE = (int)(m->iLogKE0_elec + logKE * m->iLogKE1_elec);
	bool bremsAllowed = p->ke > m->eventCutoff_phot;

	// Use random numbers and tabulated branching ratios to
	// determine which interaction takes place
	double branch = p->thread_->rng().rnd();

	if (bremsAllowed && branch < m->br10_posi[iLogKE] + logKE * m->br11_posi[iLogKE])
		DoBremsstrahlung(p, m);

	else if (branch < m->br20_posi[iLogKE] + logKE * m->br21_posi[iLogKE])
		BhabhaScatter(p, m);

	else
		AnnihilateInFlight(p);

	return e_dep;
}

// Administers a Bhabha scattering event, delegating the physics to other functions.
void mcPhysicsPositron::BhabhaScatter(mcParticle* p, const mcMediumXE* med)
{
	double cosThetaFactor = 1.0 + (TWICE_EMASS / p->ke);
	mcRng& rng = p->thread_->rng();

	// Get the energies of the positron and negatron
	double ke_posi, ke_nega;
	GetBhabhaEnergies(rng, med, p->ke, &ke_posi, &ke_nega);

	// Instantiate a new particle and assign the appropriate charges,
	// energies, and directions
	mcParticle* pNextParticle = DuplicateParticle(p);

	// Choose an arbitrary azimuthal angle
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);

	if (ke_nega < ke_posi + TWICE_EMASS)
	{
		// Negatron has less available energy and goes atop the stack
		pNextParticle->t = MCP_NEGATRON;
		pNextParticle->q = -1;
		pNextParticle->ke = ke_nega;
		double cosSquared = cosThetaFactor * ke_nega / (ke_nega + TWICE_EMASS);
		cosSquared = MIN(1.0, cosSquared);
		double cosTheta = sqrt(cosSquared);
		double sinTheta = -sqrt(1.0 - cosSquared);
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, pNextParticle->u);
		p->ke = ke_posi;
		cosSquared = cosThetaFactor * ke_posi / (ke_posi + TWICE_EMASS);
		cosSquared = MIN(1.0, cosSquared);
		cosTheta = sqrt(cosSquared);
		sinTheta = sqrt(1.0 - cosSquared);
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, p->u);
	}
	else
	{
		// Positron has less available energy and goes atop the stack
		pNextParticle->ke = ke_posi;
		double cosSquared = cosThetaFactor * ke_posi / (ke_posi + TWICE_EMASS);
		cosSquared = MIN(1.0, cosSquared);
		double cosTheta = sqrt(cosSquared);
		double sinTheta = -sqrt(1.0 - cosSquared);
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, pNextParticle->u);
		p->q = MCP_NEGATRON;
		p->q = -1;
		p->ke = ke_nega;
		cosSquared = cosThetaFactor * ke_nega / (ke_nega + TWICE_EMASS);
		cosSquared = MIN(1.0, cosSquared);
		cosTheta = sqrt(cosSquared);
		sinTheta = sqrt(1.0 - cosSquared);
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, p->u);
	}

	return;
}

void mcPhysicsPositron::GetBhabhaEnergies(mcRng& rng, const mcMediumXE* pMedium,
	double ke_in, double* ke_posi, double* ke_nega)
{
	double r; // The fraction of ke_in lost to the secondary e-

	// Compute parameters that depend only on ke_in
	double r0 = pMedium->eventCutoff_elec / ke_in;  // Minimum value of r
	double oneMinusR0 = 1.0 - r0;

	double gamma = 1.0 + ke_in / EMASS;
	double gammaSquared = SQUARE(gamma);
	double betaSquared = (gammaSquared - 1.0) / gammaSquared;

	double y = 1.0 / (gamma + 1.0);
	double ySquared = SQUARE(y);
	double x = 1.0 - 2.0 * y;
	double xSquared = SQUARE(x);

	double b4 = xSquared * x;
	double b3 = b4 + xSquared;
	double b2 = x * (3.0 + ySquared);
	double b1 = 2.0 - ySquared;

	// Loop until a value of r (between r0 and 1) is accepted
	do {
		r = r0 / (1.0 - oneMinusR0 * rng.rnd());
	} while (rng.rnd() > (1.0 - betaSquared * r*(b1 - r*(b2 - r*(b3 - r*b4)))));

	// Divide up the energy
	r = MAX(0, r);
	(*ke_nega) = r * ke_in;
	(*ke_posi) = ke_in - (*ke_nega);

	return;
}

// Administers the annihilation of a positron in flight, delegating
// the physics to other functions.
void mcPhysicsPositron::AnnihilateInFlight(mcParticle* p)
{
	double ke_in = p->ke;
	double p_in = sqrt(ke_in * (ke_in + TWICE_EMASS));
	double cosThetaFactor = p_in / ke_in;
	mcRng& rng = p->thread_->rng();

	// Get the energies of the two photons
	double e_phot1, e_phot2;
	GetAnnihilationEnergies(rng, ke_in, p_in, &e_phot1, &e_phot2);

	// Choose an arbitrary azimuthal angle
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);

	// Instantiate a new particle and assign the appropriate charges,
	// energies, and directions
	mcParticle* pNextParticle = DuplicateParticle(p);

	if (e_phot1 < e_phot2) {
		// Photon #1 has less available energy and goes atop the stack
		pNextParticle->t = MCP_PHOTON;
		pNextParticle->q = 0;
		pNextParticle->ke = e_phot1;
		double cosTheta = cosThetaFactor * (e_phot1 - EMASS) / e_phot1;
		cosTheta = MIN(1.0, cosTheta);
		double sinTheta = -sqrt((1.0 - cosTheta)*(1.0 + cosTheta));
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, pNextParticle->u);
		p->t = MCP_PHOTON;
		p->q = 0;
		p->ke = e_phot2;
		cosTheta = cosThetaFactor * (e_phot2 - EMASS) / e_phot2;
		cosTheta = MIN(1.0, cosTheta);
		sinTheta = sqrt((1.0 - cosTheta)*(1.0 + cosTheta));
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, p->u);
	}
	else {
		// Photon #2 has less available energy and goes atop the stack
		pNextParticle->t = MCP_PHOTON;
		pNextParticle->q = 0;
		pNextParticle->ke = e_phot2;
		double cosTheta = cosThetaFactor * (e_phot2 - EMASS) / e_phot2;
		cosTheta = MIN(1.0, cosTheta);
		double sinTheta = -sqrt((1.0 - cosTheta)*(1.0 + cosTheta));
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, pNextParticle->u);
		p->t = MCP_PHOTON;
		p->q = 0;
		p->ke = e_phot1;
		cosTheta = cosThetaFactor * (e_phot1 - EMASS) / e_phot1;
		cosTheta = MIN(1.0, cosTheta);
		sinTheta = sqrt((1.0 - cosTheta)*(1.0 + cosTheta));
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, p->u);
	}

	return;
}

// Computes the energy of two annihilation photons using the formulas
// on pages 269-70 of W. Heitler, The Quantum Theory of Radiation,
// 3d ed., London: Oxford University Press, 1954. Reprint: New York:
// Dover Publications, 1984. (See pages 61-66 of SLAC-265.)
void mcPhysicsPositron::GetAnnihilationEnergies(mcRng& rng, double ke_in, double p_in,
	double* e_phot1, double* e_phot2)
{
	double eAvailable = ke_in + TWICE_EMASS;
	double r; // The fraction of eAvailable going to photon #1

	// Compute parameters that depend only on ke_in and p_in
	double twiceGamma = 2.0 * (1.0 + ke_in / EMASS);
	double scaledFactor = EMASS / eAvailable;
	scaledFactor = SQUARE(scaledFactor);
	double r0 = EMASS / (eAvailable + p_in);  // Minimum value of r

	// Loop until a value of r is accepted
	do {
		r = r0 * exp(rng.rnd() * log((1.0 - r0) / r0));
	} while (rng.rnd() > 1.0 - r + scaledFactor * (twiceGamma - 1.0 / r));

	(*e_phot1) = r * eAvailable;
	(*e_phot2) = eAvailable - (*e_phot1);

	return;
}

// Creates two photons from a positron that annihilates at rest with an
// ambient electron. This function does not examine the kinetic energy
// or even the charge of the active particle: it trusts the calling function
// to ensure that the particle atop the stack is a positron that has used
// up all of its kinetic energy.
void mcPhysicsPositron::AnnihilateAtRest(mcParticle* p)
{
	mcRng& rng = p->thread_->rng();
	// Change the current particle to an appropriate photon
	p->t = MCP_PHOTON;
	p->q = 0;
	p->ke = EMASS;
	p->mfps = 0;
	GoInRandomDirection(rng.rnd(), rng.rnd(), p->u);

	// Create a second photon moving in the opposite direction:
	DuplicateParticle(p)->u *= -1.0;

	return;
}
