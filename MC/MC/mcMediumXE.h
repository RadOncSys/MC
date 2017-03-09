// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcMedium.h"

enum { BREMS_ANGLE_DEFAULT, KOCH_MOTZ };

// Pair production angular distribution options:
enum { PAIR_ANGLE_DEFAULT, LOWEST_ORDER, MOTZ_OLSEN_KOCH };

// Photoelectron angular distribution options:
enum { PHOTO_ANGLE_DEFAULT, SAUTER };

class mcMediumXE : public mcMedium
{
public:
	mcMediumXE(void);
	virtual ~mcMediumXE(void);

	void read(istream& is) override;

protected:
	static void FixStepSize(int nBins_elec, double eStep,
		double iLogKE0, double iLogKE1,
		vector<double>& stepSize0, vector<double>& stepSize1,
		vector<double>& dedx0, vector<double>& dedx1);

	void InitializeAngularDistribution();

public:
	// Photons:
	// Log energy index:
	double iLogKE0_phot;
	double iLogKE1_phot;

	// Interactions:
	double eventCutoff_phot;
	vector<double> br10_phot;         // Branching ratio 1
	vector<double> br11_phot;
	vector<double> br20_phot;         // Branching ratio 2
	vector<double> br21_phot;
	int rayleigh;                     // Switch to turn on Rayleigh scattering

	// Transport:
	vector<double> photonMFP0;        // Mean free path
	vector<double> photonMFP1;
	vector<double> raylFactor0;       // Rayleigh mean free path correction
	vector<double> raylFactor1;

	// Pair production:
	// The following data are unique to pair production:
	double aOrC[2];

	// Compton scattering:
	// No medium data needed

	// Photoelectric effect:
	double e_KEdge;

	// Rayleigh scattering:
	double  iLogKE0_rayl;
	double  iLogKE1_rayl;
	vector<double> raylQFactor0;
	vector<double> raylQFactor1;

	// Electrons (negatrons or positrons):
	// Log energy index:
	double iLogKE0_elec;
	double iLogKE1_elec;

	// Interactions:
	double eventCutoff_elec;

	// Transport:
	double radLength;          // Radiation length
	double teff0;
	double eStep;

	// Moliere multiple scattering:
	double scaledLittleB;
	double chi_cc;

	// Bremsstrahlung:
	// The following data are unique to bremsstrahlung:
	double aOrB[2];
	// The following data are shared with the pair production routines:
	double screen0[6];
	double screen1[6];
	double screen2[6];
	double screen0a[6];
	double screen1a[6];
	double screen2a[6];
	double delLimit[2];
	double del_C;     // see SLAC-265 eq. (2.7.51)
	double zFactor_angDist;  // NB: Not in the PEGS4 file

	// Negatrons:
	// Interactions:
	vector<double> br10_nega;         // Branching ratio 1
	vector<double> br11_nega;

	// Transport:
	vector<double> sigma0_nega;       // Fictitious cross section
	vector<double> sigma1_nega;
	vector<double> dedx0_nega;        // Linear energy loss rate (stopping power)
	vector<double> dedx1_nega;
	vector<double> stepSize0_nega;    // Step size
	vector<double> stepSize1_nega;

	// Moller scattering:
	double mollerThreshold;

	// Positrons:
	// Interactions:
	vector<double> br10_posi;         // Branching ratio 1
	vector<double> br11_posi;
	vector<double> br20_posi;         // Branching ratio 2
	vector<double> br21_posi;

	// Transport:
	vector<double> sigma0_posi;       // Fictitious cross section
	vector<double> sigma1_posi;
	vector<double> dedx0_posi;        // Linear energy loss rate (stopping power)
	vector<double> dedx1_posi;
	vector<double> stepSize0_posi;    // Step size
	vector<double> stepSize1_posi;

	// Bhabha scattering:
	// No medium data needed 

	// GG. Next moved from REGION to fit large XYZ requirements
	double transCutoff_phot; // Energy cutoff for photon transport
	double transCutoff_elec; // Energy cutoff for electron transport
	double ke_noRangeDisc;   // Upper kinetic energy limit for range discard 

	int bremsAngleOption;
	int pairAngleOption;
	int photoAngleOption;
};
