#include <string>
#include <vector>
#include "mcScoreTest.h"

using namespace std;

std::string typeof_(int i)
{
	{
		switch (i) {
		case 0: return "neutron";
			break;
		case 1: return "proton";
			break;
		case 2: return "deutron";
			break;
		case 3: return "triton";
			break;
		case 4: return "alpha";
			break;
		case 5: return "recoils";
			break;
		case 6: return "gamma";
			break;
		default: return "Unknown product";
			break;
		}
	}
}

void mcSpectrum::fill(double Eout, double weight, double mu, int p_type, mcRng& rng)
{
	ESpectrum.resize(3);
	int p = p_type;
	if (p_type == 6)
		p = 2;

	if (!energy_resized)
	{
		energy_resized = true;
		if (inEn / 1000000 <= 10)
			Energies.resize(100);
		else if (inEn / 1000000 <= 50)
			Energies.resize(200);
		else Energies.resize(400);
		Energies[0] = 0;
		for (int i = 1; i < Energies.size(); i++)
			Energies[i] = Energies[i - 1] + (inEn / Energies.size());
	}

	ESpectrum[p].resize(Energies.size());
	
	const double pi = 3.1415926535;
	double theta = acos(mu) * 180 / pi;
	AngSpectrum.resize(4);
	for (int i = 0; i < AngSpectrum.size(); i++)
		AngSpectrum[i].resize(2);

	if (p == 0)
	{
		if (theta < 0)
			throw exception("Theta < 0");
		else if (theta < 10)
		{
			if (Eout != 0)
			{
				AngSpectrum[0][0] += weight * Eout;
				AngSpectrum[0][1] += weight;
			}
		}
		else if (theta > 40 && theta < 50)
		{
			if (Eout != 0)
			{
				AngSpectrum[1][0] += weight * Eout;
				AngSpectrum[1][1] += weight;
			}
		}
		else if (theta > 80 && theta < 90)
		{
			if (Eout != 0)
			{
				AngSpectrum[2][0] += weight * Eout;
				AngSpectrum[2][1] += weight;
			}
		}
		else if (theta > 130 && theta < 140)
		{
			if (Eout != 0)
			{
				AngSpectrum[3][0] += weight * Eout;
				AngSpectrum[3][1] += weight;
			}
		}
	}
	/*const double pi = 3.1415926535;
	double fi = rng.rnd() * 2 * pi;
	double sr = mu * fi;*/

	for (int i = 0; i < ESpectrum[p].size(); i++)
	{
		if (Eout == 0)
			break;
		else if (Energies[i] > Eout)
		{
			ESpectrum[p][i - 1] += weight;// / Eout / 1000000) / sr);
			break;
		}
		else if (i == ESpectrum[p].size())
		{
			ESpectrum[p][i - 1] += weight;// / Eout;// / 1000000) / sr);
			break;
		}
	}
}

void mcSpectrum::dump(std::ostream& os) const
{
	os << "Incedent energy of \t" << typeof_(inP) << "\t = \t" << inEn << "\t eV" << endl;
	for (int i = 0; i < ESpectrum.size(); i++)
	{
		switch (i) {
		case 0: os << "Neutrons as a product:" << endl;
			break;
		case 1: os << "Protons as a product:" << endl;
			break;
		case 2: os << "Gamma as a product:" << endl;
			break;
		default: throw exception("Product doesn't supported for dumping spectrum");
		}
		os << "E, MeV\t" << "N" << endl;
		for (int j = 0; j < ESpectrum[i].size(); j++)
		{
			os << Energies[j] / 1000000 << "\t" << ESpectrum[i][j] << endl;
		}
		os << endl;
	}

	os << "ANGULAR DISTRIB." << endl;
	os << "0 - 10\t40 - 50\t80 - 90\t130 - 140" << endl;
	os << AngSpectrum[0][0] / AngSpectrum[0][1] << "\t" << AngSpectrum[1][0] / AngSpectrum[1][1] << "\t" << AngSpectrum[2][0] / AngSpectrum[2][1] << "\t" << AngSpectrum[3][0] / AngSpectrum[3][1] << "\t";
}
