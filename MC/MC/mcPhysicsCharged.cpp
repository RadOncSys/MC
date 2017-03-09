#include "mcPhysicsCharged.h"
#include "mcRng.h"
#include "mcDefs.h"
#include <math.h>
#include <float.h>

mcPhysicsCharged::mcPhysicsCharged(void)
{
}

mcPhysicsCharged::~mcPhysicsCharged(void)
{
}

void mcPhysicsCharged::MoliereScatter(mcRng& rng
	, double scaledLittleB
	, double chi_cc
	, geomVector3D& u
	, double eTotal
	, double betaSquared
	, double effPathLength)
{
	static const double i0BigB = 0.5714;
	static const double i1BigB = 0.21429;
	static const double bigB0[] = {
		-1.0724e+00,  3.7136e-01,  1.1396e+00,  1.4908e+00,  1.7342e+00,
		1.9233e+00,  2.0791e+00,  2.0791e+00
	};
	static const double bigB1[] = {
		2.8203e+00,  1.4560e+00,  1.1910e+00,  1.1267e+00,  1.0958e+00,
		1.0773e+00,  1.0649e+00,  1.0649e+00
	};
	static const double bigB2[] = {
		-3.5669e-01, -2.8072e-02, -5.2070e-03, -2.2565e-03, -1.2705e-03,
		-8.1806e-04, -5.7197e-04, -5.7197e-04
	};

	static const double i0G21 = 1.0;
	static const double i1G21 = 5.0;
	static const double g210[] = {
		-9.9140e-04, -9.9140e-04, -7.1017e-02, -7.3556e-02,  3.6658e-01,
		1.4498e+00,  1.4498e+00
	};
	static const double g211[] = {
		2.7672e+00,  2.7672e+00,  3.4941e+00,  3.5487e+00,  2.1162e+00,
		-5.9717e-01, -5.9717e-01
	};
	static const double g212[] = {
		-1.1544e+00, -1.1544e+00, -3.0773e+00, -3.1989e+00, -2.0311e+00,
		-3.2951e-01, -3.2951e-01
	};

	static const double i0G22 = 1.0;
	static const double i1G22 = 6.0;
	static const double g220[] = {
		-5.2593e-04, -5.2593e-04, -6.4819e-02,  3.7427e-02,  6.1955e-01,
		1.7584e+00,  2.5694e+00,  2.5694e+00
	};
	static const double g221[] = {
		1.4285e+00,  1.4285e+00,  2.2033e+00,  1.6630e+00, -6.2713e-01,
		-4.0390e+00, -6.0484e+00, -6.0484e+00
	};
	static const double g222[] = {
		-1.2670e+00, -1.2670e+00, -3.6399e+00, -2.9362e+00, -6.7859e-01,
		1.8810e+00,  3.1256e+00,  3.1256e+00
	};

	static const double i0G31 = 1.0;
	static const double i1G31 = 9.0;
	static const double g310[] = {
		4.9437e-01,  4.9437e-01,  5.3251e-01,  6.6810e-01, -3.8262e+00,
		4.2335e+00,  5.0694e+00,  1.4563e+00, -3.2852e-01, -2.2489e-01,
		-2.2489e-01
	};
	static const double g311[] = {
		1.9124e-02,  1.9124e-02, -6.1555e-01, -2.2056e+00,  2.5528e+01,
		-1.0604e+01, -1.4208e+01, -3.3275e+00,  1.2938e+00,  1.0713e+00,
		1.0713e+00
	};
	static const double g312[] = {
		1.8375e+00,  1.8375e+00,  4.5595e+00,  8.9293e+00, -3.3862e+01,
		6.6702e+00,  1.0456e+01,  2.2601e+00, -7.3254e-01, -6.1358e-01,
		-6.1358e-01
	};

	static const double i0G32 = 1.0;
	static const double i1G32 = 23.0;
	static const double g320[] = {
		2.9907e-05,  2.9907e-05,  2.5820e-03, -5.3270e-03, -6.6341e-02,
		-3.6027e-01, -2.7953e+00, -3.6091e+00,  1.2491e+01,  1.9637e+01,
		2.1692e+00, -1.6682e+01, -2.1539e+01, -1.4344e+01, -5.4990e+00,
		3.1029e+00,  6.0961e+00,  8.6179e+00,  7.5064e+00,  5.9838e+00,
		4.4959e+00,  3.2847e+00,  1.9514e+00,  4.8808e-01,  4.8808e-01
	};
	static const double g321[] = {
		4.7318e-01,  4.7318e-01,  3.5853e-01,  4.9418e-01,  1.4422e+00,
		4.7190e+00,  2.6694e+01,  3.4125e+01, -7.1103e+01, -1.1371e+02,
		-2.5019e+01,  6.2067e+01,  8.2651e+01,  5.5193e+01,  2.3874e+01,
		-4.4708e+00, -1.3670e+01, -2.0950e+01, -1.7956e+01, -1.4065e+01,
		-1.0456e+01, -7.6709e+00, -4.7505e+00, -1.6910e+00, -1.6910e+00
	};
	static const double g322[] = {
		6.5921e-01,  6.5921e-01,  1.9776e+00,  1.4528e+00, -2.2407e+00,
		-1.1380e+01, -6.0986e+01, -7.7512e+01,  9.4496e+01,  1.5794e+02,
		4.5340e+01, -5.5257e+01, -7.7065e+01, -5.0867e+01, -2.3140e+01,
		2.1318e-01,  7.2823e+00,  1.2536e+01,  1.0520e+01,  8.0342e+00,
		5.8462e+00,  4.2445e+00,  2.6452e+00,  1.0459e+00,  1.0459e+00
	};

	static const double twoMinusLogTwo = 1.306852819;
	static const double twoOverTwoMinusLogTwo = 1.530394219;

	double omega0 = scaledLittleB * effPathLength / betaSquared;

	if (omega0 > 1.0) {

		double littleB = log(omega0);

		// For the relationship between bigB and littleB, see 
		// pp. 86-87 of SLAC-265
		double bigB;
		if (littleB <= twoMinusLogTwo) {
			bigB = twoOverTwoMinusLogTwo * littleB;
		}
		else {
			int i = (int)(i0BigB + littleB * i1BigB);
			bigB = bigB0[i] + littleB * (bigB1[i] + littleB * bigB2[i]);
		}
		double inverseB = (bigB > 2.0) ? 1.0 / bigB : 0.5;

		double reducingAngle = chi_cc * sqrt(bigB * effPathLength) / (eTotal * betaSquared);

		// Compute branching ratios to be used for selecting one of the
		// three subdistribution functions: */
		double brDenominator = 1.0 + 1.75 * inverseB;
		double br1 = (1.0 - 2.0 / bigB) / brDenominator;
		double br2 = (1.0 + 0.025 * inverseB) / brDenominator;

		// Loop until an angle theta is accepted
		double rnd4;
		double theta;
		double sinTheta;
		do {
			for (;;) {
				double reducedTheta;
				for (;;) {
					int i;
					double rnd0 = rng.rnd();
					if (rnd0 <= br1) {
						// Use f0, the gaussian subdistribution
						double rnd1 = rng.rnd();
						if (rnd1 == 0.0) rnd1 = DBL_MIN;
						reducedTheta = -log(rnd1);
						reducedTheta = MAX(0, reducedTheta);
						reducedTheta = sqrt(reducedTheta);
						break;
					}
					else if (rnd0 <= br2) {
						// Use f3, the single scattering tail
						double rnd2 = rng.rnd();
						double rnd3 = rng.rnd();
						double eta = MAX(rnd2, rnd3);
						// Evaluate the rejection function g3
						i = (int)(i0G31 + eta * i1G31);
						double g31 = g310[i] + eta * (g311[i] + eta * g312[i]);
						i = (int)(i0G32 + eta * i1G32);
						double g32 = g320[i] + eta * (g321[i] + eta * g322[i]);
						double g3 = g31 + inverseB * g32;
						if (rng.rnd() <= g3) {
							reducedTheta = 1.0 / eta;
							break;
						}
					}
					else {
						// Use f2, the correction for central theta values
						reducedTheta = rng.rnd();
						// Evaluate the rejection function g2
						i = (int)(i0G21 + reducedTheta * i1G21);
						double g21 = g210[i] + reducedTheta * (g211[i] + reducedTheta * g212[i]);
						i = (int)(i0G22 + reducedTheta * i1G22);
						double g22 = g220[i] + reducedTheta * (g221[i] + reducedTheta * g222[i]);
						double g2 = g21 + inverseB * g22;
						if (rng.rnd() <= g2) break;
					}
				}
				theta = reducedTheta * reducingAngle;
				if (theta < PI) break;
			}
			sinTheta = sin(theta);
			rnd4 = rng.rnd();
		} while (SQUARE(rnd4) * theta > sinTheta);

		double cosTheta = cos(theta);
		double cosPhi, sinPhi;
		GetRandomPhi(rng.rnd(), &cosPhi, &sinPhi);
		ChangeDirection(cosTheta, sinTheta, cosPhi, sinPhi, u);
	}

	return;
}

double mcPhysicsCharged::ReducedPathLength(double stepNext, double tscat)
{
	static const double i1Path = 18.182;
	static const double path0[] = { 1.0, 1.0060, 1.0657, 1.6971, 1.6971 };
	static const double path1[] = { 0.98875, 0.78657, -0.25051, -7.56, -7.56 };
	static const double path2[] = { 2.5026, 4.2387, 8.7681, 29.946, 29.946 };

	double s = stepNext / tscat;
	int iPath = (int)(s * i1Path);
	iPath = MIN(4, iPath);
	double path = stepNext * (path0[iPath] + s * (path1[iPath] + s * path2[iPath]));

	return path;
}
