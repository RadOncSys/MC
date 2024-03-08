#include "mcPhysicsPhoton.h"
#include "mcMediumXE.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcTransport.h"
#include "mcDefs.h"

// Indices for bremsstrahlung and pair production screening functions
enum { BETHE_HEITLER, COULOMB_CORRECTED };
enum { A_BH, B_BH, C_BH, A_BHCC, B_BHCC, C_BHCC };

mcPhysicsPhoton::mcPhysicsPhoton(void)
{
}

mcPhysicsPhoton::~mcPhysicsPhoton(void)
{
}

bool mcPhysicsPhoton::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	if (p->ke <= ((const mcMediumXE&)med).transCutoff_phot || p->ke <= p->transport_->transCutoff_phot)
	{
		edep = p->ke;
		DiscardParticle(p);
		return true;
	}
	else
		return false;
}

double mcPhysicsPhoton::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	const mcMediumXE& m = (const mcMediumXE&)med;
	double logKE = log(ke);
	int iLogKE = (int)(m.iLogKE0_phot + logKE * m.iLogKE1_phot);
	if (iLogKE >= m.photonMFP0.size())
		iLogKE = m.photonMFP0.size() - 1;		//BUG??? При симуляции получил iLogKE = 200, что выходит за границы массива
	double mediumMeanFreePath = m.photonMFP0[iLogKE] + logKE * m.photonMFP1[iLogKE];

	if (m.rayleigh) {
		double rayleighFactor = m.raylFactor0[iLogKE] + logKE * m.raylFactor1[iLogKE];
		mediumMeanFreePath *= rayleighFactor;
	}

	return mediumMeanFreePath / dens;
}

double mcPhysicsPhoton::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	p->p += p->u * step;
	p->dnear -= step;
	return 0;
}

double mcPhysicsPhoton::DoInterruction(mcParticle* p, const mcMedium* med) const
{
	double e_dep = 0; // энергия, выделившаяся в результате взаимодействия (фото-эффект)

	const mcMediumXE* m = (const mcMediumXE*)med;

	double logKE = log(p->ke);
	int iLogKE = (int)(m->iLogKE0_phot + logKE * m->iLogKE1_phot);
	mcRng& rng = p->thread_->rng();

	double rayleighFactor = 1;
	if (m->rayleigh)
		rayleighFactor = m->raylFactor0[iLogKE] + logKE * m->raylFactor1[iLogKE];

	if (m->rayleigh && rng.rnd() < 1.0 - rayleighFactor)
	{
		RayleighScatter(p, m);
	}
	else
	{
		double branch = rng.rnd();

		// br1 = Pr(pair production)
		if (branch < m->br10_phot[iLogKE] + logKE * m->br11_phot[iLogKE])
			ProducePair(p, m);

		// br2 = Pr(pair production) + Pr(Compton)
		else if (branch < m->br20_phot[iLogKE] + logKE*m->br21_phot[iLogKE])
			ComptonScatter(p, m);

		else
			e_dep = DoPhotoelectric(p, m);
	}

	return e_dep;
}

void mcPhysicsPhoton::ProducePair(mcParticle* p, const mcMediumXE* pMedium)
{
	double e_in = p->ke;
	double eScaled_in = e_in / EMASS;
	double ke_elec1, ke_elec2;
	double cosTheta1, sinTheta1;
	double cosTheta2, sinTheta2, minusSinTheta2;
	mcRng& rng = p->thread_->rng();

	// Get the energies and the polar angles
	GetPairEnergies(rng, pMedium, e_in, &ke_elec1, &ke_elec2);
	GetPairAngle(rng, pMedium, eScaled_in, ke_elec1, &cosTheta1, &sinTheta1);
	GetPairAngle(rng, pMedium, eScaled_in, ke_elec2, &cosTheta2, &minusSinTheta2);
	sinTheta2 = -minusSinTheta2;

	// Choose an arbitrary azimuthal angle
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);

	// Choose randomly which particle is the negatron and which
	// is the positron
	int q_elec1, q_elec2;
	mc_particle_t t_elec1, t_elec2;
	double eAvailable_elec1, eAvailable_elec2;
	if (rng.rnd() <= 0.5) {
		t_elec1 = MCP_NEGATRON;
		t_elec2 = MCP_POSITRON;
		q_elec1 = -1;
		q_elec2 = 1;
		eAvailable_elec1 = ke_elec1;
		eAvailable_elec2 = ke_elec2 + TWICE_EMASS;
	}
	else {
		t_elec1 = MCP_POSITRON;
		t_elec2 = MCP_NEGATRON;
		q_elec1 = 1;
		q_elec2 = -1;
		eAvailable_elec1 = ke_elec1 + TWICE_EMASS;
		eAvailable_elec2 = ke_elec2;
	}

	// Instantiate a new particle and assign the appropriate charges,
	// energies, and directions
	mcParticle* pNextParticle = DuplicateParticle(p);
	if (eAvailable_elec1 < eAvailable_elec2) {
		pNextParticle->t = t_elec1;
		pNextParticle->q = q_elec1;
		pNextParticle->ke = ke_elec1;
		ChangeDirection(cosTheta1, sinTheta1, cosPhi, sinPhi, pNextParticle->u);
		p->t = t_elec2;
		p->q = q_elec2;
		p->ke = ke_elec2;
		ChangeDirection(cosTheta2, sinTheta2, cosPhi, sinPhi, p->u);
	}
	else {
		pNextParticle->t = t_elec2;
		pNextParticle->q = q_elec2;
		pNextParticle->ke = ke_elec2;
		ChangeDirection(cosTheta2, sinTheta2, cosPhi, sinPhi, pNextParticle->u);
		p->t = t_elec1;
		p->q = q_elec1;
		p->ke = ke_elec1;
		ChangeDirection(cosTheta1, sinTheta1, cosPhi, sinPhi, p->u);
	}
}

// For an incident photon energy greater than 50.0 MeV, the Coulomb-corrected
// Bethe-Heitler distribution is used; for lower incident photon energies,
// the Bethe-Heitler distribution is used. If the incident photon energy is
// less than 2.1 MeV, an approximation is made: it is assumed that one member
// of the created pair (positron or electron) receives no kinetic energy.
void mcPhysicsPhoton::GetPairEnergies(mcRng& rng, const mcMediumXE* pMedium,
	double e_in, double* ke_elec1, double* ke_elec2)
{
	double e_elec2;

	if (e_in <= 2.1)
		e_elec2 = EMASS;

	else
	{
		double r;      // Fraction of e_in that ends up with electron #2

		// We use nergies < 50 MeV
		int idist = BETHE_HEITLER;
		int ascrn = A_BH;
		int cscrn = C_BH;

		// Loop until a value of r is accepted
		for (;;)
		{
			int iscrn;

			// Decide which part of the distribution to sample from and,
			// accordingly, which screening rejection function to use
			if (rng.rnd() >= pMedium->aOrC[idist])
			{
				// Use A(delta) as the screening function
				iscrn = ascrn;

				// Sample from the subdistribution that is proportional
				// to 12 * (r-0.5)^2. By symmetry, one need only sample
				// r in the interval (0.0,0.5).
				double rnd0 = rng.rnd(),
					rnd1 = rng.rnd(),
					rnd2 = rng.rnd();
				double rndMax = MAX(rnd0, MAX(rnd1, rnd2));
				r = 0.5 * (1.0 - rndMax);
			}
			else
			{
				// Use C(delta) as the screening function
				iscrn = cscrn;
				// Sample from the subdistribution that is uniform
				r = 0.5 * rng.rnd();
			}

			if (r != 0.0)
			{
				// Make sure that A(delta) and C(delta) are positive
				double del = 1.0 / (e_in * r*(1.0 - r));
				if (del < pMedium->delLimit[idist])
				{
					double delta = pMedium->del_C * del;  // SLAC-265 eq. (2.7.52)
					double screenLevel = (delta < 1.0) ?
						pMedium->screen0[iscrn] + delta * (pMedium->screen1[iscrn] + delta *  pMedium->screen2[iscrn]) :
						pMedium->screen0a[iscrn] + pMedium->screen1a[iscrn] * log(delta + pMedium->screen2a[iscrn]);
					if (rng.rnd() <= screenLevel && r >= EMASS / e_in)
						break;
				}
			}
		}
		e_elec2 = r * e_in;
	}

	double e_elec1 = e_in - e_elec2;
	(*ke_elec1) = e_elec1 - EMASS;
	(*ke_elec2) = e_elec2 - EMASS;
}

// This function is based on NRCC report PIRS-0287R1, by Alex F. Bielajew.
//
// The three options are:
//
//     pMedia->pairAngleOption = PAIR_ANGLE_DEFAULT:
//
//         theta = 1.0 / eScaled_in;
//
//     pMedia->pairAngleOption = LOWEST_ORDER:
//
//         d(Probability)            sin(theta)
//         -------------- = -------------------------------
//            d(theta)      2*P*[E_total - P*cos(theta)]**2
//
//     pMedia->pairAngleOption = MOTZ_OLSEN_KOCH:
//
//          Uses eq. 3D-2003 of Motz, Olsen, and Koch (1969).
//
//     Note: If MOTZ_OLSEN_KOCH is chosen but the energy of the photon is
//          less than the Bethe-Heitler threshold, then the lowest-order
//          distribution is used. The default value of the threshold is
//          4.14 MeV. Users may override this with a higher value, but a
//          lower value will cause non-physical sampling.
// 
void mcPhysicsPhoton::GetPairAngle(mcRng& rng, const mcMediumXE* pMedium, double eScaled_in,
	double ke_elec, double* cosTheta, double* sinTheta)
{
	static const double betheHeitlerThreshold = 4.14 / EMASS;

	if (pMedium->pairAngleOption == LOWEST_ORDER ||
		(pMedium->pairAngleOption == MOTZ_OLSEN_KOCH && eScaled_in < betheHeitlerThreshold))
	{
		double e_elec = ke_elec + EMASS;
		double p_elec = sqrt(ke_elec * (ke_elec + TWICE_EMASS));
		double cosTemp = 1.0 - 2.0 * rng.rnd();

		(*cosTheta) = (e_elec * cosTemp + p_elec) /
			(p_elec * cosTemp + e_elec);
		double sinSquaredTemp = (1.0 - cosTemp) * (1.0 + cosTemp);
		sinSquaredTemp = MAX(0, sinSquaredTemp);
		(*sinTheta) = EMASS * sqrt(sinSquaredTemp) /
			(p_elec * cosTemp + e_elec);
	}

	else if (pMedium->pairAngleOption == MOTZ_OLSEN_KOCH) {

		double e_elec = ke_elec + EMASS;
		double eScaled_elec = e_elec / EMASS;
		double zFactor = pMedium->zFactor_angDist;

		// Calculate some parameters for repeated use
		double massFactor = 2.0 / eScaled_in;
		massFactor *= massFactor;
		double r = eScaled_elec / (eScaled_in - eScaled_elec);
		double oneOverR = 1.0 / r;
		double val1 = (1.0 + r) * (1.0 + oneOverR) / (2.0*eScaled_in);
		val1 *= val1;

		// Compute N_g (see page 5 of PIRS-0287R1)
		double y_max = eScaled_elec * PI;
		double xi_min = 1.0 / (1.0 + SQUARE(y_max));
		double g_ximin = PairRejectionFunction(xi_min, zFactor, r, oneOverR, val1);
		double xi0 = sqrt(massFactor / zFactor);
		xi0 = MAX(0.01, MAX(xi_min, MIN(0.5, xi0)));
		double alpha = 1.0 + 0.25 * log(massFactor + zFactor*SQUARE(xi0));
		double beta = 0.5 * zFactor * xi0 / (massFactor + zFactor*SQUARE(xi0));
		alpha -= beta * (xi0 - 0.5);
		double xi1 = -alpha / (3.0 * beta);
		if (alpha >= 0.0)
			xi1 += 0.5 + sqrt(0.25 + SQUARE(xi1));
		else
			xi1 += 0.5 - sqrt(0.25 + SQUARE(xi1));
		xi1 = MAX(0.01, MAX(xi_min, MIN(0.5, xi1)));
		double g_xi1 = PairRejectionFunction(xi1, zFactor, r, oneOverR, val1);
		double N_g = 1.0 / (1.02 * MAX(g_ximin, g_xi1));

		// Loop until a value of xi is accepted
		double theta, g_xi;
		do {
			double xi = rng.rnd();
			g_xi = PairRejectionFunction(xi, zFactor, r, oneOverR, val1);
			double scaledThetaSquared = (1.0 / xi) - 1.0;
			scaledThetaSquared = MAX(0, scaledThetaSquared);
			theta = sqrt(scaledThetaSquared) / eScaled_elec;
		} while (rng.rnd() > N_g * g_xi || theta > PI);

		(*cosTheta) = cos(theta);
		(*sinTheta) = sin(theta);

	}
	else {

		// Default angle selection
		double theta = 1.0 / eScaled_in;
		(*cosTheta) = cos(theta);
		(*sinTheta) = sin(theta);

	}
}

// Computes the function dG(x)/d(x) of PIRS-0287R1
double mcPhysicsPhoton::PairRejectionFunction(double x, double zFactor,
	double r, double oneOverR, double val1)
{
	double xMinusOneHalf = x - 0.5;
	double value = val1 + zFactor * x * x;
	value = 1.0 + 0.25 * log(value);
	value *= -4.0 * (r + oneOverR + 1.0 - 4.0 * xMinusOneHalf * xMinusOneHalf);
	value += 2.0 + 3.0 * (r + oneOverR);
	return value;
}

// Administers a Compton interaction, delegating the physics to other functions.
void mcPhysicsPhoton::ComptonScatter(mcParticle* p, const mcMediumXE* pMedium)
{
	//particle's parameters after interaction
	double e_phot, ke_elec;
	double cosTheta_phot, sinTheta_phot;
	double cosTheta_elec, sinTheta_elec;
	mcRng& rng = p->thread_->rng();

	// Get the energies and polar angles
	GetComptonSplit(rng, p->ke, &e_phot, &cosTheta_phot,
		&sinTheta_phot, &ke_elec, &cosTheta_elec, &sinTheta_elec);

	// Choose an arbitrary azimuthal angle
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);

	// Instantiate a new particle and assign the appropriate charges,
	// energies, and directions
	mcParticle* pNextParticle = DuplicateParticle(p);
	if (e_phot < ke_elec)
	{
		// Photon has less available energy and goes atop the stack
		pNextParticle->ke = e_phot;
		ChangeDirection(cosTheta_phot, sinTheta_phot, cosPhi, sinPhi, pNextParticle->u);
		p->t = MCP_NEGATRON;
		p->q = -1;
		p->ke = ke_elec;
		ChangeDirection(cosTheta_elec, sinTheta_elec, cosPhi, sinPhi, p->u);
	}
	else {
		// Electron has less available energy and goes atop the stack
		pNextParticle->t = MCP_NEGATRON;
		pNextParticle->q = -1;
		pNextParticle->ke = ke_elec;
		ChangeDirection(cosTheta_elec, sinTheta_elec, cosPhi, sinPhi, pNextParticle->u);
		p->ke = e_phot;
		ChangeDirection(cosTheta_phot, sinTheta_phot, cosPhi, sinPhi, p->u);
	}
}

// The calculation of the photon and electron energies and scattering
// angles have many intermediate results in common; therefore, for
// efficiency, they are all calculated together in this one function.
// Let r be the ratio of the final photon energy to the initial
// photon energy (epsilon in SLAC-265). The minimum value of r is
//     r0 = 1 / (1 + 2 * e_in / EMASS) ,
// and the maximum value of r is unity. The expression used for the
// differential cross section is proportional to
//     (1/r + r) * (1 - r * sin^2(theta)/(1+r^2)) .
// The current function samples from the first factor over the
// interval (r0, 1) and uses the second factor as a rejection function.
// 
void mcPhysicsPhoton::GetComptonSplit(mcRng& rng, double e_in,
	double* e_phot, double* cosTheta_phot, double* sinTheta_phot,
	double* ke_elec, double* cosTheta_elec, double* sinTheta_elec)
{
	// Compute parameters that depend only on e_in
	double eScaled_in = e_in / EMASS;
	double oneOverR0 = 1.0 + 2.0 * eScaled_in;
	double alpha1 = log(oneOverR0);
	double alpha2 = eScaled_in * (oneOverR0 + 1.0) / SQUARE(oneOverR0);
	double alphaSum = alpha1 + alpha2;  // Used repeatedly inside loop below
	double oneMinusR;
	double r, oneMinusCos, sinSquared;

	// Loop until an interaction is accepted
	do {
		if (alpha1 >= rng.rnd() * alphaSum)
			// Use 1/r part of distribution
			r = exp(alpha1*rng.rnd()) / oneOverR0;

		else {
			// Use linear (r) part of distribution: */
			double rPrime = rng.rnd();
			if (eScaled_in >= (eScaled_in + 1.0) * rng.rnd())
				rPrime = MAX(rPrime, rng.rnd());
			r = (1.0 + (oneOverR0 - 1.0) * rPrime) / oneOverR0;
		}
		(*e_phot) = r * e_in;
		oneMinusR = 1.0 - r;
		oneMinusCos = EMASS * oneMinusR / (*e_phot);
		sinSquared = oneMinusCos * (2.0 - oneMinusCos);
		sinSquared = MAX(0, sinSquared);
	} while (rng.rnd() > 1.0 - r * sinSquared / (1.0 + r*r));

	// Finish calculating scattering angle for photon
	(*sinTheta_phot) = sqrt(sinSquared);
	(*cosTheta_phot) = 1.0 - oneMinusCos;

	// Calculate energy and scattering angle for electron
	(*ke_elec) = e_in - (*e_phot);
	double momentumSquared = (*ke_elec) * ((*ke_elec) + TWICE_EMASS);
	if (momentumSquared <= 0.0) {
		// Avoid division by zero
		*cosTheta_elec = 0.0;
		*sinTheta_elec = -1.0;
	}
	else {
		*cosTheta_elec = ((*ke_elec) + EMASS + (*e_phot)) * oneMinusR / sqrt(momentumSquared);
		sinSquared = (1.0 - (*cosTheta_elec)) * (1.0 + (*cosTheta_elec));
		sinSquared = MAX(0, sinSquared);
		*sinTheta_elec = -sqrt(sinSquared);
	}
}

// Administers a photoelectric event, delegating the physics to other
// functions. There are two compile-time options, one to enable K-edge
// fluorescence and one to enable computation of realistic photoelectron
// angular distributions. Once enabled globally, both options can be
// turned on or off locally.
// There are two major differences between this function and other
// interaction functions: 1) this function deposits energy at the point
// of interaction; 2) the number of particles leaving the scene of the
// interaction varies: no particle if the incident photon is completely
// absorbed, one particle if only a photoelectron is generated, two
// particles if a K-edge fluorescent photon is also generated.
//
double mcPhysicsPhoton::DoPhotoelectric(mcParticle* p, const mcMediumXE* pMedium)
{
	mcRng& rng = p->thread_->rng();
	p->plast = p->p;
	// Compute the deposited energy, the photoelectron energy, and
	// the energy of the fluorescent photon (if there is one)
	double e_dep, ke_elec, e_fluor;
	GetPhotoelectricEnergies(rng.rnd(), pMedium, p->ke, &e_dep, &ke_elec, &e_fluor);

	if (ke_elec == 0) {
		// All of the incident photon energy has been absorbed
		p->ke = 0;

		// Note: The particle will be discarded at the start of the
		// next transport routine. It is not worth the miniscule speed
		// gain (if any) that would arise from discarding the particle
		// immediately here, since that would require us to reset region
		// and medium pointers inside an interaction routine.
	}
	else {
		if (e_fluor == 0) {
			// Set up the photoelectron
			p->t = MCP_NEGATRON;
			p->q = -1;
			p->ke = ke_elec;
			if (pMedium->photoAngleOption != PHOTO_ANGLE_DEFAULT && ke_elec > pMedium->transCutoff_elec)
				ChangePhotoelectronDirection(rng, ke_elec, p->u);
		}
		else {
			// Instantiate a new particle, and assign the appropriate
			// charges, energies, and directions
			mcParticle* pNextParticle = DuplicateParticle(p);

			if (e_fluor < ke_elec)
			{
				// Photon has less available energy and goes atop the stack
				pNextParticle->t = MCP_PHOTON;
				pNextParticle->q = 0;
				pNextParticle->ke = e_fluor;
				GoInRandomDirection(rng.rnd(), rng.rnd(), pNextParticle->u);
				p->t = MCP_NEGATRON;
				p->q = -1;
				p->ke = ke_elec;
				if (pMedium->photoAngleOption != PHOTO_ANGLE_DEFAULT && ke_elec > pMedium->transCutoff_elec)
					ChangePhotoelectronDirection(rng, ke_elec, p->u);
			}
			else {
				// Electron has less available energy and goes atop the stack
				pNextParticle->t = MCP_NEGATRON;
				pNextParticle->q = -1;
				pNextParticle->ke = ke_elec;
				if (pMedium->photoAngleOption != PHOTO_ANGLE_DEFAULT && ke_elec > pMedium->transCutoff_elec)
					ChangePhotoelectronDirection(rng, ke_elec, pNextParticle->u);
				p->t = MCP_PHOTON;
				p->q = 0;
				p->ke = e_fluor;
				GoInRandomDirection(rng.rnd(), rng.rnd(), p->u);
			}
		}
	}
	return e_dep;
}

// Computes the redistribution of energy in a photoelectric interaction.
// If the K-edge fluorescence option is being used, then the function
// InitializeKEdge must have been run during the initialization phase
// of the simulation. For purposes of generating K-edge fluorescent photons,
// the program treats each medium as a single element with atomic number
// Z = pMedium->z_KEdge. Note that the fluorescent photon energy is zero
// whenever the photoelectron kinetic energy is zero.
void mcPhysicsPhoton::GetPhotoelectricEnergies(double rnum, const mcMediumXE* pMedium, double e_in,
	double* e_dep, double* ke_elec, double* e_fluor)
{
	if (e_in <= pMedium->e_KEdge) {
		// Energy is below K-edge
		*e_dep = e_in;
		*ke_elec = *e_fluor = 0;
	}
	else {
		*e_fluor = 0;
		*e_dep = pMedium->e_KEdge - (*e_fluor);
		*ke_elec = e_in - pMedium->e_KEdge;
	}
}

// By default, the program assumes that a photoelectron has the same direction
// as the incident photon, but more accurate computation of the direction may
// be important in low-energy problems, e.g., Cobalt 60 gamma rays propagating
// through lead. In regions where pRegion->photoAngleOption is set to SAUTER,
// the current function computes the photoelectron direction using the Sauter
// relativistic formula. The sampling procedure is a mixed procedure that
// avoids fitting parameters to the distribution. The average acceptance
// varies from two thirds at low energy to about 0.48 at high energy. For more
// information, see: A.F. Bielajew and D.W.O. Rogers, "Photoelectron angle
// selection in the EGS4 code system," NRC Report PIRS-0058, October 1986.
//
void mcPhysicsPhoton::ChangePhotoelectronDirection(mcRng& rng, double ke_elec, geomVector3D& u)
{
	// Compute parameters that depend only on ke_elec
	double gamma = 1.0 + ke_elec / EMASS;
	double gammaSquared = SQUARE(gamma);
	double gammaFactor = (gamma + 1.0) / (2.0 * gamma);
	double beta = sqrt(ke_elec*(ke_elec + TWICE_EMASS)) / (ke_elec + EMASS);
	double alpha = 0.5 * gamma - 0.5 + 1.0 / gamma;
	double betaOverAlpha = beta / alpha;
	double xi, sinSquared, cosTheta;

	// Loop until values of cos(theta) and sin^2(theta) are accepted
	if (betaOverAlpha <= 0.2) do {
		double x = 1.0 - 2.0 * rng.rnd();
		double kappa = x + 0.5 * betaOverAlpha * (1.0 - x) * (1.0 + x);
		cosTheta = (beta + kappa) / (1.0 + beta * kappa);
		xi = 1.0 / (1.0 - beta * cosTheta);
		sinSquared = (1.0 - cosTheta) * (1.0 + cosTheta);
		sinSquared = MAX(0, sinSquared);
	} while (rng.rnd() > gammaFactor * xi * sinSquared);
	else do {
		double x = 1.0 - 2.0 * rng.rnd();
		xi = gammaSquared * (1.0 + alpha *
			(sqrt(1.0 + betaOverAlpha * (2.0*x + betaOverAlpha)) - 1.0));
		cosTheta = (1.0 - 1.0 / xi) / beta;
		sinSquared = (1.0 - cosTheta) * (1.0 + cosTheta);
		sinSquared = MAX(0, sinSquared);
	} while (rng.rnd() > gammaFactor * xi * sinSquared);

	double sinTheta = sqrt(sinSquared);
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);
	ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, u);
}

void mcPhysicsPhoton::GetRayleighAngle(mcRng& rng, const mcMediumXE* pMedium,
	double e_phot, double* cosTheta, double* sinTheta)
{
	static const double rayleighConstant = EMASS_SQUARED / SQUARE(20.60744);
	// Compute a parameter that depends only on e_phot
	double eFactor = 0.5 / (e_phot * e_phot);
	double cosSquared;

	for (;;) {
		double rnd0 = rng.rnd();
		int iLogKE = (int)(pMedium->iLogKE0_rayl + rnd0 * pMedium->iLogKE1_rayl);
		double qFactor = pMedium->raylQFactor0[iLogKE] +
			rnd0 * pMedium->raylQFactor1[iLogKE];
		double qSquared = qFactor * rayleighConstant;
		*cosTheta = 1.0 - qSquared * eFactor;
		if (ABS(*cosTheta) > 1.0) continue;
		cosSquared = SQUARE(*cosTheta);
		if (rng.rnd() <= 0.5 * (1.0 + cosSquared)) break;
	}

	double sinSquared = 1.0 - cosSquared;
	sinSquared = MAX(0, sinSquared);
	*sinTheta = sqrt(sinSquared);
}

// Administers a Rayleigh interaction, delegating the physics to
// other functions.
void mcPhysicsPhoton::RayleighScatter(mcParticle* p, const mcMediumXE* pMedium)
{
	mcRng& rng = p->thread_->rng();
	// Get the polar angle
	double cosTheta, sinTheta;
	GetRayleighAngle(rng, pMedium, p->ke, &cosTheta, &sinTheta);

	// Choose an arbitrary azimuthal angle
	double cosPhi, sinPhi;
	GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);

	// Set the new direction
	ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, p->u);
	p->plast = p->p;
}
