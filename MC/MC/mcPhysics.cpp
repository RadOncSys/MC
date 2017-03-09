#include "mcPhysics.h"
#include "mcDefs.h"
#include "mcParticle.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

mcPhysics::mcPhysics(void)
{
}

mcPhysics::~mcPhysics(void)
{
}

// Returns the cosine and sine of a random angle between zero and two pi.
void mcPhysics::GetRandomPhi(double rnum, double* cosPhi, double* sinPhi)
{
	double phi = rnum * TWOPI;
	*cosPhi = cos(phi);
	*sinPhi = sin(phi);
}

// Generates new direction cosines from the scattering angles and old direc-
// tion cosines. Theta is the polar scattering angle in the transport frame,
// and phi is the azimuthal scattering angle in the transport frame. In the
// laboratory frame, let psi be the polar angle and delta the azimuthal angle
// of the direction in which the particle is initially travelling. The trans-
// formation from the laboratory frame to the transport frame is performed by
// first rotating an angle delta about the z-axis and then rotating an angle
// psi about the y-axis. In the present case, the scattering angles in the
// transport frame are known, and it is necessary to transform the correspond-
// ing direction cosines back into the laboratory frame. This reverse trans-
// formation is performed by first rotating an angle -psi about the y-axis
// and then rotating an angle -delta about the z-axis. The final direction
// cosines in the laboratory frame are thus obtained by means of the matrix
// multiplication:
//
//  /                        \  /                    \  /                 \
//  | cosDelta  -sinDelta  0 |  |  cosPsi  0  sinPsi |  | sinTheta cosPhi |
//  |                        |  |                    |  |                 |
//  | sinDelta   cosDelta  0 |  |    0     1    0    |  | sinTheta sinPhi |
//  |                        |  |                    |  |                 |
//  |    0          0      1 |  | -sinPsi  0  cosPsi |  | cosTheta        |
//  \                        /  \                    /  \                 /
//
void mcPhysics::ChangeDirection(double cosTheta, double sinTheta,
	double cosPhi, double sinPhi, geomVector3D& u)
{
	static const double minSinPsiSquared = 1.0e-10;
	double sinPsiCosDelta = u.x();
	double sinPsiSinDelta = u.y();
	double cosPsi = u.z();
	double sinThetaCosPhi = sinTheta * cosPhi;
	double sinThetaSinPhi = sinTheta * sinPhi;
	double sinPsiSquared = sinPsiCosDelta * sinPsiCosDelta + sinPsiSinDelta * sinPsiSinDelta;

	if (sinPsiSquared < minSinPsiSquared) {
		if (cosPsi > 0)
			u.set(sinThetaCosPhi, sinThetaSinPhi, cosTheta);
		else
			u.set(-sinThetaCosPhi, -sinThetaSinPhi, -cosTheta);
	}
	else {
		double sinPsi = sqrt(sinPsiSquared);
		double cosDelta = sinPsiCosDelta / sinPsi;
		double sinDelta = sinPsiSinDelta / sinPsi;
		u.set((cosPsi * cosDelta * sinThetaCosPhi) - (sinDelta * sinThetaSinPhi) + (sinPsiCosDelta * cosTheta),
			(cosPsi * sinDelta * sinThetaCosPhi) + (cosDelta * sinThetaSinPhi) + (sinPsiSinDelta * cosTheta),
			-(sinPsi * sinThetaCosPhi) + (cosPsi * cosTheta));
	}
}

// Computes a set of direction cosines from an isotropic distribution.
void mcPhysics::GoInRandomDirection(double rnum1, double rnum2, geomVector3D& u)
{
	double cosTheta = 1.0 - 2.0 * rnum1;
	double sinSquared = (1.0 - cosTheta) * (1.0 + cosTheta);
	sinSquared = MAX(0, sinSquared);
	double sinTheta = sqrt(sinSquared);
	double cosPhi, sinPhi;
	GetRandomPhi(rnum2, &cosPhi, &sinPhi);
	u.set(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
}

mcParticle* mcPhysics::DuplicateParticle(mcParticle* p)
{
	p->plast = p->p;
	return p->thread_->DuplicateParticle();
}

void mcPhysics::DiscardParticle(mcParticle* p)
{
	p->thread_->RemoveParticle();
}
