#include "mcPhysicsElectron.h"
#include "mcMediumXE.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcTransport.h"
#include "mcDefs.h"
#include <float.h>

enum { BETHE_HEITLER, COULOMB_CORRECTED };
enum { A_BH, B_BH, C_BH, A_BHCC, B_BHCC, C_BHCC };

mcPhysicsElectron::mcPhysicsElectron(void)
{
}

mcPhysicsElectron::~mcPhysicsElectron(void)
{
}

bool mcPhysicsElectron::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	if (p->ke <= ((const mcMediumXE&)med).transCutoff_elec || p->ke <= p->transport_->transCutoff_elec)
	{
		edep = p->ke;
		DiscardParticle(p);
		return true;
	}
	else
		return false;
}

double mcPhysicsElectron::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	const mcMediumXE& m = (const mcMediumXE&)med;
	double logKE = log(ke);
	int iLogKE = int(m.iLogKE0_elec + logKE * m.iLogKE1_elec);
	double sigma = (m.sigma0_nega[iLogKE] + logKE * m.sigma1_nega[iLogKE])* dens;
	return (sigma > 0.0) ? 1 / sigma : DBL_MAX;
}

double mcPhysicsElectron::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	static const double minBetaSquared = 1.0e-8;
	static const double twoOverEsubS = 0.094315;

	const mcMediumXE& m = (const mcMediumXE&)med;
	double e_dep;

	// «десь нет проблемы пересечени€ границы. ќна уже решена раньше.
	// ¬се, что требуетс€ - это определить шаг, на который нам следует переместить частицу
	// согласно выбранной стратегии моделировани€ непрерывных потерь и рассе€ни€.
	// «атем, перемещаем на рассто€ние, равное минимумму из запрошенного и шага непреывных потерь.
	// ¬ конце коректируем направление скорости согласно теории рассе€ни€ и обновл€ем переменную шага
	// (что нужно дл€ моделировани€ дискретных событий).

	double logKE = log(p->ke);
	int iLogKE = (int)(m.iLogKE0_elec + logKE * m.iLogKE1_elec);
	double dedx = p->regDensityRatio * (m.dedx0_nega[iLogKE] + logKE * m.dedx1_nega[iLogKE]);

	// —итуаци€, когда частица заведомо не выйдет из области.
	// Ўаг не может быть больше 
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

	// –еальный шаг выбираетс€ как минимальный из шага до дискретного событи€, 
	// шага непрерывных потерь и рассто€ние, на котором вли€ние рассе€ни€ не станвитс€ критичным дл€ модели.
	double stepSize = (m.stepSize0_nega[iLogKE] + logKE * m.stepSize1_nega[iLogKE]) / p->regDensityRatio;
	step = MIN(stepSize, MIN(step, 0.3 * tscat));

	// ƒлина трека больше рассто€ни€ между точками и св€зана рассеивающей поспособностью.
	double pathLength = ReducedPathLength(step, tscat);

	// ѕотери в треке расчитываютс€ в два шага, сначала дл€ dedx текущей энергии, затем средней.
	e_dep = dedx * pathLength;
	e_dep = MIN(e_dep, p->ke);

	logKE = log(p->ke - 0.5 * e_dep);
	iLogKE = int(m.iLogKE0_elec + logKE * m.iLogKE1_elec);
	if (iLogKE < 0) iLogKE = 0;  // вблизи энергии поглощени€ возможна проблема индексов
	dedx = p->regDensityRatio * (m.dedx0_nega[iLogKE] + logKE * m.dedx1_nega[iLogKE]);
	e_dep = dedx * pathLength;
	e_dep = MIN(e_dep, p->ke);
	p->p += p->u * step;

	// Change direction in accordance with the Moliere theory
	MoliereScatter(p->thread_->rng(), m.scaledLittleB, m.chi_cc, p->u,
		p->ke + EMASS, betaSquared,
		pathLength * p->regDensityRatio);

	p->plast = p->p;
	p->mfps = 0;
	p->ke -= e_dep;
	p->dnear -= step;
	return e_dep;
}

double mcPhysicsElectron::DoInterruction(mcParticle* p, const mcMedium* med) const
{
	const mcMediumXE* m = (const mcMediumXE*)med;
	double e_dep = 0;

	bool bremsAllowed = p->ke > m->eventCutoff_phot;
	bool mollerAllowed = p->ke > m->mollerThreshold;
	mcRng& rng = p->thread_->rng();

	if (!bremsAllowed)
		MollerScatter(p, m);
	else
	{
		if (!mollerAllowed)
			DoBremsstrahlung(p, m);
		else
		{
			double logKE = log(p->ke);
			int iLogKE = (int)(m->iLogKE0_elec + logKE * m->iLogKE1_elec);
			double br1 = m->br10_nega[iLogKE] + logKE * m->br11_nega[iLogKE];
			if (rng.rnd() < br1)
				DoBremsstrahlung(p, m);
			else
				MollerScatter(p, m);
		}
	}

	return e_dep;
}

// Administers a Moller scattering event, delegating the physics to
// other functions. Inasmuch as the kinetic energy transfer is never
// more than half the incident kinetic energy, the kinetic energy of
// the incident electron must be greater than or equal to twice
// eventCutoff_elec. 

void mcPhysicsElectron::MollerScatter(mcParticle* p, const mcMediumXE* med)
{
	double cosThetaFactor = 1.0 + (TWICE_EMASS / p->ke);
	mcRng& rng = p->thread_->rng();

	// Get the energies of the two electrons
	double ke_elec1, ke_elec2;
	GetMollerEnergies(rng, med, p->ke, &ke_elec1, &ke_elec2);

	// Instantiate a new particle and assign the appropriate charges, energies, and directions
	mcParticle* pNextParticle = DuplicateParticle(p);

	// The function GetMollerEnergies() assures that ke_elec2 is less than
	// or equal to ke_elec1, so the electron with energy ke_elec2 goes atop the stack
	pNextParticle->ke = ke_elec2;
	double cosSquared = cosThetaFactor * ke_elec2 / (ke_elec2 + TWICE_EMASS);
	cosSquared = MIN(1.0, cosSquared);
	double cosTheta = sqrt(cosSquared);
	double sinTheta = -sqrt(1.0 - cosSquared);
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);
	ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, pNextParticle->u);
	p->ke = ke_elec1;
	cosSquared = cosThetaFactor * ke_elec1 / (ke_elec1 + TWICE_EMASS);
	cosSquared = MIN(1.0, cosSquared);
	cosTheta = sqrt(cosSquared);
	sinTheta = sqrt(1.0 - cosSquared);
	ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, p->u);
	p->plast = p->p;

	return;
}

// Administers a bremsstrahlung event, delegating the physics to
// other functions.
// NOTE: There is another version that allows bremsstrahlung
// splitting. Be sure to correct any bugs in both versions.

void mcPhysicsElectron::DoBremsstrahlung(mcParticle* p, const mcMediumXE* med)
{
	// Get the electron and photon energies
	double ke_in = p->ke;
	double e_phot = BremsPhotonEnergy(p->thread_->rng(), med, ke_in);
	double ke_elec = ke_in - e_phot;

	// Instantiate a new particle and assign the appropriate charges,
	// energies, and directions
	mcParticle* pNextParticle = DuplicateParticle(p);
	double eAvailable_elec = (pNextParticle->q == 1) ? ke_elec + TWICE_EMASS : ke_elec;
	if (eAvailable_elec < e_phot) {
		// Electron has less available energy and goes atop the stack
		pNextParticle->ke = ke_elec;
		SetUpBremsPhoton(p, med, ke_in, e_phot);
	}
	else {
		// Photon has less available energy and goes atop the stack
		SetUpBremsPhoton(pNextParticle, med, ke_in, e_phot);
		p->ke = ke_elec;
	}

	return;
}

// Computes the energies of the two electrons leaving a Moller scattering
// event. By definition, ke_elec2 is the smaller of the two energies.

void mcPhysicsElectron::GetMollerEnergies(mcRng& rng, const mcMediumXE* med,
	double ke_in, double* ke_elec1, double* ke_elec2)
{
	double keScaled_in = ke_in / EMASS;
	double eAboveThreshold = ke_in - med->mollerThreshold;
	double r;       // The fraction of ke_in lost to the secondary e-
	double ratio;

	// Compute parameters that depend only on ke_in
	double gamma = 1.0 + keScaled_in;   //  e_in/EMASS  (special relativity)
	double gammaSquared = SQUARE(gamma);

	double g2 = SQUARE(keScaled_in) / gammaSquared;
	double g3 = (2.0 * keScaled_in + 1.0) / gammaSquared;
	double gmax = (1.0 + 1.25 * g2);

	// Loop until a value of r is accepted
	do {
		r = med->eventCutoff_elec /
			(ke_in - eAboveThreshold * rng.rnd());
		ratio = r / (1.0 - r);
	} while (gmax * rng.rnd() > (1.0 + g2*SQUARE(r) + ratio * (ratio - g3)));

	// Divide up the energy
	(*ke_elec2) = r * ke_in;
	(*ke_elec1) = ke_in - (*ke_elec2);

	return;
}

// For electron total energies greater than 50.0 MeV, the Coulomb-corrected
// Bethe-Heitler distribution is used; otherwise, the Bethe-Heitler distri-
// bution is used.
double mcPhysicsElectron::BremsPhotonEnergy(mcRng& rng, const mcMediumXE* pMedium, double ke_in)
{
	static const double oneOverLog2 = 1.442695041;
	static const double oneOverTwiceLog2 = 0.721347520;
	double  e_in = ke_in + EMASS;
	double  e_phot;

	int idist;
	int ascrn;
	int bscrn;
	if (e_in <= 50.0) {
		idist = BETHE_HEITLER;
		ascrn = A_BH;
		bscrn = B_BH;
	}
	else {
		idist = COULOMB_CORRECTED;
		ascrn = A_BHCC;
		bscrn = B_BHCC;
	}

	// Compute number of subdistributions needed to produce photons of minimum
	// energy, in case the (1-r)/r part of the distribution is used
	int n_brems = (int)(oneOverLog2 * log(e_in / pMedium->eventCutoff_phot));

	// NOTE: The preceding equation agrees with the EGS4 code, but it agrees
	// with eq. (2.7.70) of SLAC-265 only if the right-hand side is increased
	// by unity. */

	// Loop until a photon energy is accepted
	for (;;) {

		// Decide which part of the distribution to sample from and,
		// accordingly, which screening rejection function to use
		// (see the "if" conditions in eqs. (2.7.91-92) of SLAC 265)

		int iscrn;
		double r; // ratio of photon energy to incident electron total energy

		if (rng.rnd() * (n_brems * pMedium->aOrB[idist] + 0.5) > 0.5)
		{

			// Use A(delta) as the screening function
			iscrn = ascrn;

			// Sample from the (1-r)/r part of the distribution,
			// choosing a subdistribution at random from the set
			// (0,1,2,...,NBREMS-1) (eqs. (2.7.76-77) of SLAC-265)
			double p = pow(0.5, double(int(n_brems * rng.rnd())));

			// All subdistributions are sampled by first sampling from:
			//    1.0/log(2.0)                   if 0.0 <= r <  0.5
			//   (1.0/log(2.0)) * (1.0-r)/r      if 0.5 <= r <= 1.0
			// and then taking r = r * p. In the initial sampling, the
			// probability that r is less than 0.5 is 1/(2*log(2))
			if (rng.rnd() < oneOverTwiceLog2) {
				// Sample from the r < 0.5 part
				r = 0.5 * rng.rnd();
			}
			else {
				// Sample from the r >= 0.5 part
				double rnd0;
				do {
					rnd0 = rng.rnd();
					double rnd1 = rng.rnd();
					double rnd2 = rng.rnd();
					double rndMax = MAX(rnd1, rnd2);
					r = 1.0 - 0.5 * rndMax;
				} while (rnd0 > 0.5 / r);
			}
			r *= p;

		}
		else {

			// Use B(delta) as the screening function
			iscrn = bscrn;

			// Sample from the 2*r part of the distribution
			r = MAX(rng.rnd(), rng.rnd());

		}

		// Calculate the energy of the bremsstrahlung photon: */
		e_phot = e_in * r;
		if (e_phot >= pMedium->eventCutoff_phot && e_phot <= ke_in) {
			/* Make sure that A(delta) and B(delta) are positive: */
			double del = r / (e_in - e_phot);
			if (del < pMedium->delLimit[idist]) {
				double delta = pMedium->del_C * del;  /* SLAC-265 eq. (2.7.52) */
				double screenLevel;
				if (delta < 1.0) {
					screenLevel = pMedium->screen0[iscrn] + delta *
						(pMedium->screen1[iscrn] +
							delta * pMedium->screen2[iscrn]);
				}
				else {
					screenLevel = pMedium->screen0a[iscrn] +
						pMedium->screen1a[iscrn] *
						log(delta + pMedium->screen2a[iscrn]);
				}
				if (rng.rnd() <= screenLevel) break;
			}
		}
	}

	return e_phot;
}

// This is function g(x) from eq. (4) of PIRS-0203. The parameters val1,
// val2, and val3 are r-dependent quantities which have been computed
// in advance to save time.
double mcPhysicsElectron::BremsRejectionFunction(double x, double zFactor, double r, double val1, double val2, double val3)
{
	double temp = x + 1.0;
	temp *= temp;
	double g = val2 + (4.0*r * x / temp - val1) * (4.0 + log(val3 + zFactor / temp));

	// NOTE: The above expression for g(x) agrees with EGS4 and appears to
	// be correct in that it avoids negative return values; however, the
	// sign of val2 is opposite to that found in eq. (4) of PIRS-0203.

	return g;
}

// This function is based on NRCC Report PIRS-0203, by Alex F. Bielajew,
// Rahde Mohan, and Chen-Shou Chui.
// The two options are:
//     pMedia->bremsAngleOption = BREMS_ANGLE_DEFAULT:
//         theta = EMASS / e_in;
//     pMedia->bremsAngleOption = KOCH_MOTZ:
//          Uses eq. 2BS of Koch and Motz (1959).

double mcPhysicsElectron::BremsPhotonAngle(mcRng& rng, const mcMediumXE* pMedium, double ke_in, double e_phot)
{
	double e_in = ke_in + EMASS;
	double theta;  // scattering angle 

	if (pMedium->bremsAngleOption == KOCH_MOTZ) {

		double gamma = e_in / EMASS;
		double zFactor = pMedium->zFactor_angDist;

		double y_max = gamma * PI;
		double x_max = SQUARE(y_max);

		// Calculate some parameters for repeated use by the rejection function g(x)
		double r = 1.0 - e_phot / e_in;
		double val1 = 1.0 + SQUARE(r);
		double val2 = 3.0 * val1 - 2.0 * r;
		double val3 = (1.0 - r) / (2.0 * r * gamma);
		val3 *= val3;

		// Compute N_r (see page 5 of PIRS-0203)
		double g_xmin = BremsRejectionFunction(0.0, zFactor, r, val1, val2, val3);
		double g_xmid = BremsRejectionFunction(1.0, zFactor, r, val1, val2, val3);
		double g_xmax = BremsRejectionFunction(x_max, zFactor, r, val1, val2, val3);
		// Normalization factor (see eq. (6) of PIRS-0203)
		double N_r = 1.0 / MAX(g_xmin, MAX(g_xmid, g_xmax));

		// Loop until a value of x is accepted
		double x;
		double g_x;    // Value of rejection function g(x)
		do {
			// x equals the value under the radical in equation (5) of PIRS-0203
			double rnd0 = rng.rnd();
			x = rnd0 / (1.0 - rnd0 + 1.0 / x_max);
			// Compute g(x): */
			g_x = BremsRejectionFunction(x, zFactor, r, val1, val2, val3);
		} while (rng.rnd() > N_r * g_x);

		// Convert the accepted value of x into an angle in radians
		theta = sqrt(x) / gamma;
	}
	else {

		// Default angle selection
		theta = EMASS / e_in;
	}

	return theta;
}

void mcPhysicsElectron::SetUpBremsPhoton(mcParticle* pPhoton, const mcMediumXE* med, double ke_in, double e_phot)
{
	// Assign the correct charge and energy
	pPhoton->t = MCP_PHOTON;
	pPhoton->q = 0;
	pPhoton->ke = e_phot;

	// Establish the new direction
	mcRng& rng = pPhoton->thread_->rng();
	double theta = BremsPhotonAngle(rng, med, ke_in, e_phot);
	double cosTheta = cos(theta);
	double sinTheta = sin(theta);
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);
	ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, pPhoton->u);

	return;
}
