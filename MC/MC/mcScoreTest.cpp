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

void mcSpectrum::fill(double Eout, double weight, int p_type)
{
	ESpectrum.resize(3);
	int p = p_type;
	if (p_type == 6)
		p = 2;
	if (inEn / 1000000 > 100)
		ESpectrum[p].resize(150);
	else ESpectrum[p].resize(round(inEn / 1000000));
	Eout /= 1000000;
	for (int i = 0; i < ESpectrum[p].size(); i++)
	{
		if (Eout == 0)
			break;
		else if (i > Eout)
		{
			ESpectrum[p][i - 1] += weight;
			break;
		}
		else if (i == ESpectrum[p].size())
		{
			ESpectrum[p][i - 1] += weight;
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
		os << "E, MeV\t" << "E * n" << endl;
		for (int j = 0; j < ESpectrum[i].size(); j++)
		{
			os << j << "\t" << ESpectrum[i][j] << endl;
		}
		os << endl;
	}
}
