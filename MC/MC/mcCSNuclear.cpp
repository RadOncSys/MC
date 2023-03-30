#include "mcCSNuclear.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>

using namespace std;

void mcCSNuclearForAngleSpectrum::Clear()
{
	SourceBinUps.clear();
	SourceSpectrum.clear();
}

void mcCSNuclearParticleEnergy::Clear()
{
	Energy = 0;
	TotalCrossSection = 0;
	ProtonCrossSection = 0;
	NeutronCrossSection = 0;

	ProtonAngles.clear();
	NeutronAngles.clear();
}

void mcCSNuclear::Load(const char* fname, const char* ename)
{
	// Separators
	static const string energySeparator("+++++++++++++++++++++++++++++++");
	static const string evalueSeparator("Cross section summary information for");
	static const string tcsvalueSeparator("Nonelastic cross section =");
	static const string protonsSeparator("*protons*");
	static const string neutronsSeparator("*neutrons*");
	static const string protonCrossSectionSeparator("proton production cross section (mb) =");
	static const string neutronCrossSectionSeparator("neutron production cross section (mb) =");

	ifstream isIcru(fname);
	if (isIcru.fail())
		throw exception((string("Can't open Proton data file: ") + ename).c_str());

	// База данных изотопа
	ElementName = ename;

	// Сечения для энергии падающей частицы
	mcCSNuclearParticleEnergy csForEnergy;
	mcCSNuclearForAngleSpectrum particleAngle;
	mcCSNuclearForAngleSpectrum neutronAngle;

	// Читаем строки текста одну за другой и выбираем нужную информацию
	string line, s1, s2, s3, s4;
	std::getline(isIcru, line, '\n');

	// Состояния указыват в каком месте парсинга мы находимся и потому как интерпитируем строки
	bool isNewEnergyPrepare = false;
	bool isProtons = false;
	bool isNeutrons = false;
	int angleCount = 0;

	while (!isIcru.fail())
	{
		// Начало новой энергии
		if (line.find(energySeparator) != string::npos)
		{
			csForEnergy.Clear();
			isNewEnergyPrepare = true;
		}
		else if (isNewEnergyPrepare)
		{
			if (line.find(evalueSeparator) != string::npos)
			{
				csForEnergy.Energy = atof(&line[evalueSeparator.size()]);
			}
			else if (line.find(tcsvalueSeparator) != string::npos)
			{
				csForEnergy.TotalCrossSection = atof(&line[tcsvalueSeparator.size()]);
				isNewEnergyPrepare = false;
			}
		}

		// proton spectra
		else if (line.find(protonsSeparator) != string::npos)
		{
			csForEnergy.ProtonAngles.clear();
			isProtons = true;
			angleCount = 0;
		}
		else if (isProtons)
		{
			if (line.find("ANGLE(deg)") != string::npos)
			{
				std::vector<std::string> ss;
				GetStringArray(line, ss, " ");
				angleCount = (int)ss.size() - 2;
				csForEnergy.ProtonAngles.resize(angleCount);
				for (int i = 0; i < angleCount; i++)
					csForEnergy.ProtonAngles[i].Angle = atof(ss[i + 1].c_str());
			}
			else if (line.find("ENERGY") != string::npos || line.find("------") != string::npos)
			{
			}
			else if (line.find(protonCrossSectionSeparator) != string::npos)
			{
				vector<string> ss;
				GetStringArray(line, ss, "=");
				csForEnergy.ProtonCrossSection = atof(ss[1].c_str());
				isProtons = false;
			}
			else
			{
				vector<string> ss;
				GetStringArray(line, ss, " ");
				for (int i = 0; i < (int)csForEnergy.ProtonAngles.size(); i++)
				{
					csForEnergy.ProtonAngles[i].SourceBinUps.push_back(atof(ss[1].c_str()));
					csForEnergy.ProtonAngles[i].SourceSpectrum.push_back(atof(ss[i + 2].c_str()));
				}
			}

		}

		// neutron spectra
		else if (line.find(neutronsSeparator) != string::npos)
		{
			csForEnergy.NeutronAngles.clear();
			isNeutrons = true;
			angleCount = 0;
		}
		else if (isNeutrons)
		{
			if (line.find("ANGLE(deg)") != string::npos)
			{
				vector<string> ss;
				GetStringArray(line, ss, " ");
				angleCount = (int)ss.size() - 2;
				csForEnergy.NeutronAngles.resize(angleCount);
				for (int i = 0; i < angleCount; i++)
					csForEnergy.NeutronAngles[i].Angle = atof(ss[i + 1].c_str());
			}
			else if (line.find("ENERGY") != string::npos || line.find("------") != string::npos)
			{
			}
			else if (line.find(neutronCrossSectionSeparator) != string::npos)
			{
				vector<string> ss;
				GetStringArray(line, ss, "=");
				csForEnergy.NeutronCrossSection = atof(ss[1].c_str());
				isNeutrons = false;
			}
			else
			{
				vector<string> ss;
				GetStringArray(line, ss, " ");
				for (int i = 0; i < (int)csForEnergy.NeutronAngles.size(); i++)
				{
					csForEnergy.NeutronAngles[i].SourceBinUps.push_back(atof(ss[1].c_str()));
					csForEnergy.NeutronAngles[i].SourceSpectrum.push_back(atof(ss[i + 2].c_str()));
				}
			}
		}

		else if (line.find("**************************************************************************************************************************") != string::npos)
		{
			Energies.push_back(csForEnergy);
		}
		getline(isIcru, line, '\n');
	}
}

void mcCSNuclear::Clear()
{
	Energies.clear();
}
