#include "mcEndfP.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>
#include <iostream>

using namespace std;

double mcEndfRecord::ParseValue(const char* s, int n)
{
	double f = atof(s);
	int i = 0;
	for (; i < n; i++)
		if (s[i] == '+' || s[i] == '-') break;
	if (i < n)
		f *= pow(10.0, atof(s + i));
	return f;
}

void mcEndfCrossSectionTable::Load(istream& is)
{
	string line, s1, s2, s3, s4;
	mcEndfRecord record;
	int pointCount = 0;

	getline(is, line, '\n');

	while (!is.fail())
	{
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		// Сечения 
		if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '2')
		{
			ninterpolations = atoi(record.c[4]);
			if (ninterpolations > 1)
				throw exception("UUnexpected multiple interpolation types");
		}
		else if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '3')
		{
			npoints = atoi(record.c[0]);
			interpolationType = atoi(record.c[1]);
			Energies.resize(npoints, 0);
			Values.resize(npoints, 0);
		}
		else
		{
			for (int ii = 0; ii < 6; ii += 2)
			{
				if (pointCount < npoints)
				{
					Energies[pointCount] = mcEndfRecord::ParseValue(record.c[ii], 11);
					Values[pointCount] = mcEndfRecord::ParseValue(record.c[ii + 1], 11);
					pointCount++;
				}
			}
			if (pointCount == npoints)
				break; // конец таблицы, прерываем while
		}

		getline(is, line, '\n');
	}
}

void mcEndfEANuclearCrossSectionTable::mLoad(istream& is)
{
	string line, s1, s2, s3, s4;
	mcEndfRecord record;
	int pointCount = 0;

	getline(is, line, '\n');

	while (!is.fail())
	{
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		// Мультиплетности 
		if (record.c[2][10] == ' ')
		{
			n_energypoints = atoi(record.c[0]);
			interpolationType = atoi(record.c[1]);
			Energies.resize(n_energypoints, 0);
			Multiplicities.resize(n_energypoints, 0);
		}
		else
		{
			for (int ii = 0; ii < 6; ii += 2)
			{
				if (pointCount < n_energypoints)
				{
					Energies[pointCount] = mcEndfRecord::ParseValue(record.c[ii], 11);
					Multiplicities[pointCount] = mcEndfRecord::ParseValue(record.c[ii + 1], 11);
					pointCount++;
				}
			}
		}
		if (pointCount == n_energypoints)
			break;
		getline(is, line, '\n');
	}
}

void mcEndfEANuclearCrossSectionTable::Load(std::istream& is)
{
	string line;
	mcEndfRecord record;
	int pointCount = 0, i = 0;
	bool is_ea_read = false;
	getline(is, line, '\n');
	::memcpy(&record, line.c_str(), 80);
	//чтение энерго-угловых параметров
	while (!is.fail())
	{
		
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		if (!is_ea_read)
		{
			n_energypoints = atoi(record.c[0]);
			ninterpolations = atoi(record.c[1]);
			EA_par.resize(n_energypoints);
			getline(is, line, '\n');
			is_ea_read = true;
		}
		else
		{
			for (i = 0; i < n_energypoints; i++)
			{
				npoints_ang = atoi(record.c[5]);
				EA_par[i].resize(npoints_ang);
				
				for (int k = 0; k < EA_par[i].size(); k++)
					EA_par[i][k].resize(3);
				
				for (int j = 0; j < npoints_ang; j++)
				{												
					if (j % 2 == 0)
					{
						getline(is, line, '\n');
						::memcpy(&record, line.c_str(), 80);
					}
					EA_par[i][j][0] = mcEndfRecord::ParseValue(record.c[(j % 2)*3], 11);
					EA_par[i][j][1] = mcEndfRecord::ParseValue(record.c[(j % 2) * 3 + 1], 11);
					EA_par[i][j][2] = mcEndfRecord::ParseValue(record.c[(j % 2) * 3 + 2], 11);
				}
				if (pointCount < n_energypoints - 1)
				{
					getline(is, line, '\n');
					::memcpy(&record, line.c_str(), 80);
				}
				pointCount++;
			}
		}
		if (pointCount == n_energypoints)
			break;
	}
}

mcEndfProduct::mcEndfProduct()
{
}

mcEndfProduct::~mcEndfProduct()
{
	for (int i = 0; i < EANuclearCrossSections.size(); i++)
	{
		delete EANuclearCrossSections[i];
	}
}

void mcEndfProduct::Load(std::istream& is)
{
	// Читаем строки текста одну за другой и выбираем нужную информацию
	string line, s1, s2, s3, s4;
	int i = 0;
	getline(is, line, '\n');

	// Состояния указыват в каком месте парсинга мы находимся и потому как интерпитируем строки
	bool isInData = false;
	bool is_multiplicity_read = false;
	int energy_count = 0;
	int pointCount = 0;

	mcEndfRecord record;

	while (!is.fail())
	{
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		if (!is_multiplicity_read)
		{
			ZAP = atoi(record.c[0]);
			AWP = atof(record.c[1]);

			switch (ZAP) {
			
				case 1:
					product_type = particle_type::neutron;
					break;
				case 1001:
					product_type = particle_type::proton;
					break;
				case 1002:
					product_type = particle_type::deutron;
					break;
				case 1003:
					product_type = particle_type::triton;
					break;
				case 2004:
					product_type = particle_type::alpha;
					break;
				case 0:
					product_type = particle_type::gamma;
					break;
				default:
					product_type = particle_type::recoils;
					break;
			}
			auto multiplicities = new mcEndfEANuclearCrossSectionTable();
			multiplicities->mLoad(is);
			NuclearMultiplicity.push_back(multiplicities);
			is_multiplicity_read = true;
		}
		else
		{
			auto energyangle = new mcEndfEANuclearCrossSectionTable();
			energyangle->Load(is);
			EANuclearCrossSections.push_back(energyangle);
			break;																	//Разобраться с брейками!! Считывать NK в начале, чтобы считывать все продукты
		}

 		getline(is, line, '\n');
	}

}

void mcEndfCrossSectionTable::dump(std::ostream& os) const
{
	os << "NPoints = \t" << npoints << endl;
	os << "InterpolationType = \t" << interpolationType << endl;
	os << endl;
	os << "Energy\tValue" << endl;

	for (int i = 0; i < Energies.size(); i++)
		os << Energies[i] << "\t" << Values[i] << endl;
}

// dump y_i


mcEndfP::mcEndfP()
{
}

mcEndfP::~mcEndfP()
{
	for (int i = 0; i < Products.size(); i++)
	{
		delete Products[i];
	}
}

void mcEndfP::Load(const char* fname, const char* ename)
{
	// Separators
	static const string beginSeparator("*** C O N T E N T S ***");
	
	ifstream isEndf(fname);
	if (isEndf.fail())
		throw exception((string("Can't open Proton data file: ") + ename).c_str());
	ElementName = ename;

	// Читаем строки текста одну за другой и выбираем нужную информацию
	string line, s1, s2, s3, s4;
	getline(isEndf, line, '\n');

	// Состояния указыват в каком месте парсинга мы находимся и потому как интерпитируем строки
	bool isInData = false;
	int pointCount = 0;

	mcEndfRecord record;

	while (!isEndf.fail())
	{
		if(line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		// Начало новой энергии
		if (!isInData)
		{
			if (line.find(beginSeparator) != string::npos)
				isInData = true;
		}

		// Последняя строка файла. Прерываем не дожидаясь ошибки.
		else if (record.Stblt[0] == '-' && record.Stblt[1] == '1')
			break;

		// Длина строки из девяток является разделителем между таблицами
		else if (string(record.LineNumber, 5) == "99999")
		{
			pointCount = 0;
		}

		// Используем только MF=3 (сечения реакций) / MT=5 (сумма всех реакций за исключением отдельно оговоренных)
		// и     MF=6 (энерго-угловые распределени) / MT=5

		// Сечения суммы эластичных рассеяний и ядерных реакций
		else if (record.MF[0] == ' ' && record.MF[1] == '3' && 
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '2')
		{
			TotalCrossSections.Load(isEndf);
		}

		// Сечения ядерных реакций
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')
		{
			NuclearCrossSections.Load(isEndf);
		}

		// Энерго-угловые распределения
		else if (record.MF[0] == ' ' && record.MF[1] == '6' && 
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')
		{
			if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
			{
				auto product = new mcEndfProduct();
				product->Load(isEndf);
				Products.push_back(product);
			}
		}

		getline(isEndf, line, '\n');
	}
}

void mcEndfP::Clear()
{
	//Energies.clear();
}

void mcEndfP::dumpTotalCrossections(ostream& os) const
{
	os << endl;
	os << "Dump total proton crossections for element = \t" << ElementName << endl;
	os << "---------------------------------------------------------------" << endl;
	os << endl;
	TotalCrossSections.dump(os);

	os << endl;
	os << "Dump nuclear proton crossections for element = \t" << ElementName << endl;
	os << "--------------------------------------------------------------" << endl;
	os << endl;
	NuclearCrossSections.dump(os);

	os << endl;
	os << "Dump EA proton crossections for element = \t" << ElementName << endl;
	os << "---------------------------------------------------------------" << endl;
	os << endl;
	//EANuclearCrossSections.dump(os);
}


