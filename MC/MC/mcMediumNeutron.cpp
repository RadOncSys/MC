#include "mcMediumNeutron.h"
#include "../geometry/text.h"
#include <iostream>
#include "mcPhysicsCommon.h"

#ifndef InverseRadiationLengthNeutron
#define InverseRadiationLengthNeutron InverseRadiationLength_DahlApproximationNeutron
#endif InverseRadiationLengthNeutron

double InverseRadiationLength_DahlApproximationNeutron(const double A, const double Z)
{
	return Z * (Z + 1) * log(287 / sqrt(Z)) / (716.4 * A);
}

double InverseRadiationLengthNeutron(const double* A, const double* Z, const double* w, const int n)
{
	double s = 0;
	for (int i = 0; i < n; i++)
		s += InverseRadiationLengthNeutron(A[i], Z[i]) * w[i];
	return s;
}

double rmsNuclearRadiusNeutron(int At)
{
	double ac =
		(At == 1) ? 0.84 :
		(At == 2) ? 2.71 :
		(At == 3) ? 1.78 :
		(At == 4) ? 1.63 :
		((At >= 6) && (At <= 14)) ? 2.4 :
		(0.82 * pow(At, 1.0 / 3.0) + 0.58);	// At>=16

	// Предполагаем, что в формуле 4.84 ошибка - лишний корень
	return sqrt(SQUARE(ac) - 0.64);
}

double sigmaTripathiLightNeutron(int Ap, int Zp, int At, int Zt, double KE) {
	double r0 = 100.0 * 1.1E-15; // в сантиметрах, чтобы сразу возвращать рез. в см2
	double Z3 = 1.0 / 3.0;
	double Ap3 = pow(Ap, Z3);
	double At3 = pow(At, Z3);
	double Apt3 = Ap3 + At3;
	double Ecm = ECenterOfMass(At * PMASS, Ap * PMASS, Ap * PMASS + KE); // приближённое вычисление масс
	double Ecm3 = pow(Ecm, Z3);
	double SL = 1.2 + 1.6 * (1 - exp(-KE / 15.0));
	double X1 = 2.83 - 3.1E-2 * At + 1.7E-4 * SQUARE(At); // для протонов
	double Xm = 1.0 - X1 * exp(-KE / (X1 * SL));
	double D = 1.85 + 0.16 / (1.0 + exp((500.0 - KE) / 200.0));
	double T1 = 18; // для протонов, хотя есть подозрение, что 18 
	double CE = D * (1.0 - exp(-KE / T1)) - 0.292 * exp(-KE / 792.0) * cos(0.229 * pow(KE, 0.453));// было ошибочно 0.291,
	double S = Ap3 * At3 / Apt3;
	double deltaE = (1.85 + 0.16 / Ecm3) * S - CE + 0.91 * (At - 2 * Zt) * Zp / (At * Ap);
	double rTrms = rmsNuclearRadiusNeutron(At); // root mean square radius in fm 
	double rPrms = rmsNuclearRadiusNeutron(Ap); // аппроксимация Willson et al. 1991 [VK add 06.08]
	double rT = 1.29 * rTrms;
	double rP = 1.29 * rPrms;
	double R = rP + rT + 1.2 * (Apt3) / Ecm3;
	double B = 1.44 * Zp * Zt / R;
	double RC =
		(At == 1 && Zt == 1) ? 7.0 : // p+p - моё предположение - линейная экстарполяция. Если это всё вообще применимо для pp
		(At == 2 && Zt == 1) ? 13.5 : // p+d
		(At == 3 && Zt == 2) ? 21.0 : // p+3He
		(At == 4 && Zt == 2) ? 27.0 : // p+4He
		(At == 6 && Zt == 3) ? 2.2 : // p+4He
		1.0; // по аналогии с другими формулами где нет этого множителя и общей тенденцией
	double fff = r0 * (Apt3 + deltaE);
	double sigmaTL = PI * SQUARE(fff) * (1 - RC * B / Ecm) * Xm;
	sigmaTL = (sigmaTL > 0.0) ? sigmaTL : 0.0;
	return sigmaTL;
}

void coeff_calcNeutron(const vector<double>& s, vector<double>& a, vector<double>& b)
{
	int n = (int)s.size();
	a.assign(n, 0.0); b.assign(n, 0.0);
	for (int i = 0; i < n - 1; i++)
	{
		a[i] = s[i + 1] - s[i];
		b[i] = s[i] - a[i] * (i + 1);
	}
	a[n - 1] = a[n - 2];
	b[n - 1] = b[n - 2];
	return;
}

mcMediumNeutron::mcMediumNeutron()
	:transCutoff_neutron(0.1)
{
}

mcMediumNeutron::~mcMediumNeutron(void)
{
}

void mcMediumNeutron::read(istream& is)
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
	int i, ne = atoi(ParseLine(line, "NE").c_str());
	this->elements_.resize(ne);
	for (i = 0; i < ne; i++) {
		// Line
		getline(is, line, '\n');
		if (is.fail()) return;

		mcElement& e = this->elements_[i];
		strcpy_s(e.atomicSymbol, 3, ParseLine(line, "ASYM").c_str());
		e.atomicNumber = atoi(ParseLine(line, "Z").c_str());
		e.atomicMass = atof(ParseLine(line, "A").c_str());
		e.partsByNumber = atof(ParseLine(line, "PZ").c_str());
	}

	for (;;) 
	{ 
		getline(is, line, '\n');
		if (is.fail()) return;
		if (line.substr(0, 5) == "TABLE")break;
	}

	getline(is, line, '\n');
	if (is.fail()) return;
	vector<string> ss; int nl;
	nl = GetStringArray(line, ss, "\t"); // разбираем заголовок на столбцы
	int ikE = -1, idEdx = -1;
	for (int i = 0; i < nl; i++) { 
		if (ss[i] == "kE") { ikE = i; }
		if (ss[i] == "dE/dx") { idEdx = i; }
	}

	if ((ikE == -1) || (idEdx == -1)) { throw std::exception("kE or dE/dx data not found"); }
	getline(is, line, '\n'); 
	vector<double>kE, dEdx;
	getline(is, line, '\n');
	i = 0;
	while (line.substr(0, 3) != "END") {
		vector<double> data;
		int nd = GetFloatArray(line, data);
		if (nd != nl) { throw std::exception("Wrong Table"); }
		kE.push_back(data[ikE]); dEdx.push_back(data[idEdx] * density_);
		getline(is, line, '\n');
		i++;
	}
	for (i = 0; i < kE.size(); i++) {
		if (kE[i] != double(i + 1)) { throw std::exception("kE is not [1,2,...]"); }
	};
	coeff_calcNeutron(dEdx, dedx1_proto, dedx0_proto);

	gdEdxStragglingGaussVarianceConstPart();
	gRadiationLength();
	gSigmaInelastic();

	status_ = LOADED;
}

const double mcMediumNeutron::AtomicWeight() const
{
	double A = 0.0;
	for (vector<mcElement>::const_iterator el = elements_.begin(); el != elements_.end(); el++)
		A += el->atomicMass * el->partsByNumber;
	return A;
}

const double mcMediumNeutron::gdEdxStragglingGaussVarianceConstPart()
{
	dEdxStragglingGaussVarianceConstPart_ = 0.0;
	for (vector<mcElement>::iterator el = elements_.begin(); el != elements_.end(); el++)
		dEdxStragglingGaussVarianceConstPart_ += el->partsByNumber * el->atomicNumber;
	dEdxStragglingGaussVarianceConstPart_ *= 0.3 * SQUARE(EMASS) * density_ / AtomicWeight();
	return dEdxStragglingGaussVarianceConstPart_;
};

const double mcMediumNeutron::gRadiationLength()
{
	radLength = 0.0;
	for (vector<mcElement>::iterator el = elements_.begin(); el != elements_.end(); el++) {
		radLength += InverseRadiationLengthNeutron(el->atomicMass, el->atomicNumber) * el->partsByNumber * el->atomicMass;
	}
	radLength = AtomicWeight() / (radLength * density_);
	return radLength;
}

const void	mcMediumNeutron::gSigmaInelastic(int Ap, int Zp)
{
	double S;
	vector<double>sigma_in;

	for (int i = 0; i < kEmax(); i++) 
	{
		S = 0.0; // длина свободного пробега
		for (vector<mcElement>::iterator el = elements_.begin(); el != elements_.end(); el++) 
		{
			S += sigmaTripathiLightNeutron(Ap, Zp,
				ROUND(el->atomicMass), ROUND(el->atomicNumber),
				i + 1) * el->partsByNumber;
		}
		sigma_in.push_back(S * NAVOGADRO * density_ / AtomicWeight());
	}
	coeff_calcNeutron(sigma_in, sigma1_proto, sigma0_proto);
}
