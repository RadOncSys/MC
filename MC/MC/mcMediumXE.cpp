#include "mcMediumXE.h"
#include "mcDefs.h"
#include "../geometry/text.h"
#include <math.h>

mcMediumXE::mcMediumXE(void)
	: iLogKE0_phot(0), iLogKE1_phot(0), eventCutoff_phot(0)
	, rayleigh(0), e_KEdge(0), iLogKE0_rayl(0), iLogKE1_rayl(0), iLogKE0_elec(0), iLogKE1_elec(0)
	, eventCutoff_elec(0), radLength(0), teff0(0), eStep(0.2)
	, scaledLittleB(0), chi_cc(0), del_C(0)
	, zFactor_angDist(0), mollerThreshold(0)
	, transCutoff_phot(0.05)
	, transCutoff_elec(0.2)
	, ke_noRangeDisc(0)
	, bremsAngleOption(KOCH_MOTZ)
	, pairAngleOption(MOTZ_OLSEN_KOCH)
	, photoAngleOption(SAUTER)
{
}

mcMediumXE::~mcMediumXE(void)
{
}

void mcMediumXE::read(istream& is)
{
	const double distanceUnit = 1.0;

	status_ = FAILED;

	string line, s1, s2;
	vector<double> a;
	vector<int> ia;

	// Line
	getline(is, line, '\n');
	if (is.fail()) return;

	// Тип среды, плотность
	string eletype;
	GetTwoStringsFromLine(line, eletype, s2);
	line = s2;
	this->density_ = atof(ParseLine(line, "RHO").c_str());

	// Элементарный состав
	int i, j, ne = atoi(ParseLine(line, "NE").c_str());
	this->elements_.resize(ne);
	for (i = 0; i < ne; i++)
	{
		// Line
		getline(is, line, '\n');
		if (is.fail()) return;

		mcElement& e = this->elements_[i];
		strcpy_s(e.atomicSymbol, 3, ParseLine(line, "ASYM").c_str());
		e.atomicNumber = atoi(ParseLine(line, "Z").c_str());
		e.atomicMass = atof(ParseLine(line, "A").c_str());
		e.partsByNumber = atof(ParseLine(line, "PZ").c_str());
	}

	// Line
	getline(is, line, '\n');
	if (is.fail()) return;
	if (GetFloatArray(line, a) != 5) throw std::exception("Five parameters expected in the line");

	this->radLength = a[0];
	this->eventCutoff_elec = a[1];
	this->eventCutoff_phot = a[2];
	// Максимальные энергии фотонов и электронов, для которых имеются сечения.
	// В оригинальной версии использовались для контроля удовлетворения PEGS файла задаче.
	double eMax_elec = a[3];
	double eMax_phot = a[4];

	// Note: The following lines adjust for the fact that PEGS reports total
	// energies whereas this program carries kinetic energies. */
	this->eventCutoff_elec -= EMASS;
	this->mollerThreshold = 2.0 * this->eventCutoff_elec;

	// Line
	getline(is, line, '\n');
	if (is.fail()) return;
	if (GetIntArray(line, ia) != 9)
		throw std::exception("Nine parameters expected in the line");

	int nBins_phot = ia[1];
	int nBins_elec = ia[3];
	this->rayleigh = ia[7];

	// Read bremsstrahlung data:
	// 8 строк объединяем в одну, что будет соответсвовать всему блоку
	line.clear();
	for (i = 0; i < 8; i++) {
		// Line
		getline(is, s1, '\n');
		if (is.fail()) return;
		line += s1;
	}
	if (GetFloatArray(line, a) != 36) throw std::exception("36 parameters expected");

	j = 0;
	for (i = 0; i < 6; i++)
	{
		this->screen0[i] = a[j++];
		this->screen1[i] = a[j++];
		this->screen2[i] = a[j++];
		this->screen0a[i] = a[j++];
		this->screen1a[i] = a[j++];
		this->screen2a[i] = a[j++];
	}

	// Line
	getline(is, line, '\n');
	if (is.fail()) return;
	getline(is, s1, '\n');
	if (is.fail()) return;
	line += s1;
	if (GetFloatArray(line, a) != 7) throw std::exception("7 parameters expected");

	this->del_C = a[0];
	this->aOrB[0] = a[1];
	this->aOrC[0] = a[2];
	this->delLimit[0] = a[3];
	this->aOrB[1] = a[4];
	this->aOrC[1] = a[5];
	this->delLimit[1] = a[6];

	// Read electron data:
	// Line
	getline(is, line, '\n');
	if (is.fail()) return;
	getline(is, s1, '\n');
	if (is.fail()) return;
	line += s1;
	if (GetFloatArray(line, a) != 6) throw std::exception("6 parameters expected");

	this->teff0 = a[1];
	this->scaledLittleB = a[2];
	this->chi_cc = a[3];
	this->iLogKE0_elec = a[4];
	this->iLogKE1_elec = a[5];

	this->sigma0_nega.resize(nBins_elec);
	this->sigma1_nega.resize(nBins_elec);
	this->sigma0_posi.resize(nBins_elec);
	this->sigma1_posi.resize(nBins_elec);
	this->dedx0_nega.resize(nBins_elec);
	this->dedx1_nega.resize(nBins_elec);
	this->dedx0_posi.resize(nBins_elec);
	this->dedx1_posi.resize(nBins_elec);
	this->br10_nega.resize(nBins_elec);
	this->br11_nega.resize(nBins_elec);
	this->br10_posi.resize(nBins_elec);
	this->br11_posi.resize(nBins_elec);
	this->br20_posi.resize(nBins_elec);
	this->br21_posi.resize(nBins_elec);
	this->stepSize0_nega.resize(nBins_elec);
	this->stepSize1_nega.resize(nBins_elec);
	this->stepSize0_posi.resize(nBins_elec);
	this->stepSize1_posi.resize(nBins_elec);

	int n = (nBins_elec * 16 + 4) / 5;
	line.clear();
	for (i = 0; i < n; i++) {
		// Line
		getline(is, s1, '\n');
		if (is.fail()) return;
		line += s1;
	}
	if (GetFloatArray(line, a) != nBins_elec * 16)
		throw std::exception("nBins_elec * 16 parameters expected");

	j = 0;
	for (i = 0; i < nBins_elec; i++)
	{
		this->sigma0_nega[i] = a[j++];
		this->sigma1_nega[i] = a[j++];
		this->sigma0_posi[i] = a[j++];
		this->sigma1_posi[i] = a[j++];
		this->dedx0_nega[i] = a[j++];
		this->dedx1_nega[i] = a[j++];
		this->dedx0_posi[i] = a[j++];
		this->dedx1_posi[i] = a[j++];
		this->br10_nega[i] = a[j++];
		this->br11_nega[i] = a[j++];
		this->br10_posi[i] = a[j++];
		this->br11_posi[i] = a[j++];
		this->br20_posi[i] = a[j++];
		this->br21_posi[i] = a[j++];
		this->stepSize0_nega[i] = a[j++];
		this->stepSize1_nega[i] = a[j++];
		this->stepSize0_posi[i] = this->stepSize0_nega[i];
		this->stepSize1_posi[i] = this->stepSize1_nega[i];
	}

	// Read photon data:
	// Line
	getline(is, line, '\n');
	if (is.fail()) return;
	if (GetFloatArray(line, a) != 3) throw std::exception("3 parameters expected");

	this->e_KEdge = a[0];
	this->iLogKE0_phot = a[1];
	this->iLogKE1_phot = a[2];

	this->photonMFP0.resize(nBins_phot);
	this->photonMFP1.resize(nBins_phot);
	this->br10_phot.resize(nBins_phot);
	this->br11_phot.resize(nBins_phot);
	this->br20_phot.resize(nBins_phot);
	this->br21_phot.resize(nBins_phot);

	n = (nBins_phot * 6 + 4) / 5;
	line.clear();
	for (i = 0; i < n; i++) {
		// Line
		getline(is, s1, '\n');
		if (is.fail()) return;
		line += s1;
	}
	if (GetFloatArray(line, a) != nBins_phot * 6)
		throw std::exception("nBins_phot * 6 parameters expected");

	j = 0;
	for (i = 0; i < nBins_phot; i++)
	{
		this->photonMFP0[i] = a[j++];
		this->photonMFP1[i] = a[j++];
		this->br10_phot[i] = a[j++];
		this->br11_phot[i] = a[j++];
		this->br20_phot[i] = a[j++];
		this->br21_phot[i] = a[j++];
	}

	if (this->rayleigh > 0)
	{
		this->raylFactor0.resize(nBins_phot);
		this->raylFactor1.resize(nBins_phot);

		// Line
		getline(is, line, '\n');
		if (is.fail()) return;
		int nBins_rayl = atoi(line.c_str());

		// Line
		getline(is, line, '\n');
		if (is.fail()) return;
		if (GetFloatArray(line, a) != 2) throw std::exception("2 parameters expected");
		this->iLogKE0_rayl = a[0];
		this->iLogKE1_rayl = a[1];

		this->raylQFactor0.resize(nBins_rayl);
		this->raylQFactor1.resize(nBins_rayl);

		// QFactors
		n = (nBins_rayl * 2 + 4) / 5;
		line.clear();
		for (i = 0; i < n; i++) {
			// Line
			getline(is, s1, '\n');
			if (is.fail()) return;
			line += s1;
		}
		if (GetFloatArray(line, a) != nBins_rayl * 2)
			throw std::exception("nBins_rayl * 2 parameters expected");

		j = 0;
		for (i = 0; i < nBins_rayl; i++)
		{
			this->raylQFactor0[i] = a[j++];
			this->raylQFactor1[i] = a[j++];
		}

		// Factors
		n = (nBins_phot * 2 + 4) / 5;
		line.clear();
		for (i = 0; i < n; i++) {
			// Line
			getline(is, s1, '\n');
			if (is.fail()) return;
			line += s1;
		}
		if (GetFloatArray(line, a) != nBins_phot * 2)
			throw std::exception("nBins_phot * 2 parameters expected");

		j = 0;
		for (i = 0; i < nBins_phot; i++)
		{
			this->raylFactor0[i] = a[j++];
			this->raylFactor1[i] = a[j++];
		}
	}

	// Scale parameters delivered by PEGS in units of radiation length:
	this->radLength /= distanceUnit;
	double distanceFactor = this->radLength;
	double inverseDistanceFactor = 1.0 / distanceFactor;
	for (i = 0; i < nBins_phot; i++) {
		this->photonMFP0[i] *= distanceFactor;
		this->photonMFP1[i] *= distanceFactor;
	}
	for (i = 0; i < nBins_elec; i++) {
		this->sigma0_nega[i] *= inverseDistanceFactor;
		this->sigma1_nega[i] *= inverseDistanceFactor;
		this->sigma0_posi[i] *= inverseDistanceFactor;
		this->sigma1_posi[i] *= inverseDistanceFactor;
		this->dedx0_nega[i] *= inverseDistanceFactor;
		this->dedx1_nega[i] *= inverseDistanceFactor;
		this->dedx0_posi[i] *= inverseDistanceFactor;
		this->dedx1_posi[i] *= inverseDistanceFactor;
		this->stepSize0_nega[i] *= distanceFactor;
		this->stepSize1_nega[i] *= distanceFactor;
		this->stepSize0_posi[i] *= distanceFactor;
		this->stepSize1_posi[i] *= distanceFactor;
	}
	this->teff0 *= distanceFactor;
	this->scaledLittleB *= inverseDistanceFactor;
	this->chi_cc *= sqrt(inverseDistanceFactor);

	// Correct for the fact that PEGS assumes matrix indices
	// begin with one, whereas they begin with zero in C:
	this->iLogKE0_elec -= 1.0;
	this->iLogKE0_phot -= 1.0;
	if (this->rayleigh > 0)
		this->iLogKE0_rayl -= 1.0;

	// Fix step sizes (note that the differing dedx arrays for negatrons
	// and positrons here give rise to differing stepSize arrays):
	FixStepSize(nBins_elec, this->eStep,
		this->iLogKE0_elec, this->iLogKE1_elec,
		this->stepSize0_nega, this->stepSize1_nega,
		this->dedx0_nega, this->dedx1_nega);

	FixStepSize(nBins_elec, this->eStep,
		this->iLogKE0_elec, this->iLogKE1_elec,
		this->stepSize0_posi, this->stepSize1_posi,
		this->dedx0_posi, this->dedx1_posi);

	// Initialize bremsstrahlung and pair production angular distribution:
	InitializeAngularDistribution();

	status_ = LOADED;
}

// Changes the electron step size arrays so that each step results in a
// fixed energy loss. For a detailed discussion, see "Low Energy Electron
// Transport with EGS," Nucl. Inst. Meth. A, 227, 535--48 (1984). 
void mcMediumXE::FixStepSize(int nBins_elec, double eStep,
	double iLogKE0, double iLogKE1,
	vector<double>& stepSize0, vector<double>& stepSize1,
	vector<double>& dedx0, vector<double>& dedx1)
{
	double logKE = -iLogKE0 / iLogKE1;
	double dedx = dedx0[0] + logKE * dedx1[0];
	double saveStep = eStep * exp(logKE) / dedx;

	for (int i = 0; i < nBins_elec - 1; i++) {

		logKE += 1.0 / iLogKE1;   // Fixed interval between log KE values
		dedx = dedx0[i + 1] + logKE * dedx1[i + 1];
		double step = eStep * exp(logKE) / dedx;

		// Solve the equations:
		//     saveStep = stepSize0 + savelogKE * stepSize1 ,
		//         step = stepSize0 +     logKE * stepSize1
		stepSize1[i] = (step - saveStep) * iLogKE1;
		stepSize0[i] = step - logKE * stepSize1[i];

		saveStep = step;
	}

	// Last table entry applies only to last energy:
	stepSize0[nBins_elec - 1] = stepSize0[nBins_elec - 2];
	stepSize1[nBins_elec - 1] = stepSize1[nBins_elec - 2];
}

// Initializes angular distribution data for bremsstrahlung or pair
// production. The quantity zFactor_angDist equals [(1/111)*z_eff^(1/3)]^2,
// where z_eff is defined in equation (7) of NRCC report PIRS-0287 by Alex F. Bielajew.
void mcMediumXE::InitializeAngularDistribution()
{
	static const double oneOver111Squared = 1.0 / (111.0 * 111.0);
	static const double oneThird = 1.0 / 3.0;

	double norm = this->zFactor_angDist = 0.0;

	for (int i = 0; i < (int)this->elements_.size(); i++) {
		this->zFactor_angDist += this->elements_[i].partsByNumber *
			this->elements_[i].atomicNumber *
			(1.0 + this->elements_[i].atomicNumber);
		norm += this->elements_[i].partsByNumber;
	}
	this->zFactor_angDist /= norm;
	this->zFactor_angDist = oneOver111Squared * pow(this->zFactor_angDist, oneThird);
}
