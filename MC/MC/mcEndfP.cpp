#include "mcEndfP.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>
#include <iostream>

using namespace std;

double mcEndfRecord::ParseValue(const char* s, int n)
{
	std::string s1 = s;
	s1.erase(n);
	double f = atof(s);
	int i = 0;
	for (; i < n; i++)											
		if (s[i] == '+' || s[i] == '-' && i != 0) break;  //Исправлен парсинг upd. 17.05.23 by GK
	if (i < n)
	{
		s1.erase(0, i);
		f *= pow(10.0, std::stoi(s1));
	}
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
				//	cout << record.c[ii + 1] << endl;
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
	bool is_e_points_read = false;

	getline(is, line, '\n');

	while (!is.fail())
	{
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		// Мультиплетности 
		if (!is_e_points_read)
		{
			n_energypoints = atoi(record.c[0]);
			interpolationType = atoi(record.c[1]);
			Energies.resize(n_energypoints, 0);
			Multiplicities.resize(n_energypoints, 0);
			is_e_points_read = true;
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
	LANG = atoi(record.c[2]);

	//чтение энерго-угловых параметров
	getline(is, line, '\n');

	if (LANG == 2)
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
					NA = atoi(record.c[3]);
					if (NA == 2)
					{
						throw exception((string("NA = 2 for Kalbach-Mann doesn't supported ") + line).c_str());
					}
					npoints_out = mcEndfRecord::iStrCrop(record.c[5], 11);
					EA_par[i].resize(npoints_out);

					for (int k = 0; k < EA_par[i].size(); k++)
						EA_par[i][k].resize(3);

					for (int j = 0; j < npoints_out; j++)
					{
						if (j % 2 == 0)
						{
							getline(is, line, '\n');
							::memcpy(&record, line.c_str(), 80);
						}
						EA_par[i][j][0] = mcEndfRecord::ParseValue(record.c[(j % 2) * 3], 11);
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
	else if (LANG == 1)
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
					int c = 0, c1 = 0; //counters
					NA = atoi(record.c[3]);
					npoints_out = mcEndfRecord::iStrCrop(record.c[5], 11);
					EA_par[i].resize(npoints_out);

					for (int k = 0; k < EA_par[i].size(); k++)
						EA_par[i][k].resize(NA + 2);

					for (int j = 0; j < npoints_out; j++)
					{
						if (j == 0)
						{
							getline(is, line, '\n');
							::memcpy(&record, line.c_str(), 80);
						}
						EA_par[i][j][0] = mcEndfRecord::ParseValue(record.c[c1], 11);
						c = c1 + 1;
						for (int ii = 0; ii < NA + 1; ii++)
						{
							EA_par[i][j][ii + 1] = mcEndfRecord::ParseValue(record.c[c], 11);
							if (c % 5 == 0 && j!=npoints_out - 1)
							{
								getline(is, line, '\n');
								::memcpy(&record, line.c_str(), 80);
								c = 0;
							}
							else
								c++;
						}
						c1 = c;
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
	else
		throw exception("This LANG type isn't supported by this code");
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

std::string mcEndfProduct::typeof(int i)
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

	auto energyangle = new mcEndfEANuclearCrossSectionTable();

	mcEndfRecord record;

	while (!is.fail())
	{
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		if (!is_multiplicity_read)
		{
			//long double test = atof(record.c[0]);
			ZAP = round(mcEndfRecord::ParseValue(record.c[0], 11));
			AWP = mcEndfRecord::ParseValue(record.c[1], 11);

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
			energyangle->mLoad(is);
			is_multiplicity_read = true;
		}
		else
		{
			energyangle->Load(is);
			EANuclearCrossSections.push_back(energyangle);
			break;																
		}
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

void mcEndfEANuclearCrossSectionTable::dump(std::ostream& os) const
{
	os << "Multiplicity table:" << endl;
	os << "# \t" << "Energy \t" << "Multiplicity" << endl;
	for (int i = 0; i < n_energypoints; i++)
	{
		os << i << "\t" << Energies[i] << "\t" << Multiplicities[i] << endl;
	}
	if (LANG == 2)
	{
		os << "---------------------KALBACH-MAHN REPRESENTATION------------------" << endl;
		for (int i = 0; i < n_energypoints; i++)
		{
			os << "#" << i << "\t" << "Incedent energy = \t" << Energies[i] << endl;
			os << "Out energy \t " << "f_0 \t " << "r \t" << endl;
			for (int j = 0; j < EA_par[i].size(); j++)
			{
				os << EA_par[i][j][0] << "\t" << EA_par[i][j][1] << "\t" << EA_par[i][j][2] << endl;
			}
		}
		os << endl;
	}
	else if (LANG == 1)
	{
		os << "---------------------LEGANDRE REPRESENTATION------------------" << endl;
		for (int i = 0; i < n_energypoints; i++)
		{
			os << "#" << i << "\t" << "Incedent energy = \t" << Energies[i] << endl;
			for (int j = 0; j < EA_par[i].size(); j++)
			{
				os << "Out energy \t ";
				for (int k = 0; k < EA_par[i][j].size() - 1; k++)
				{
					os << "f_" << k << "\t";
				}
				os << endl;
				for (int k = 0; k < EA_par[i][j].size(); k++)
				{
					os << EA_par[i][j][k] << "\t";
				}
				os << endl;
			}
		}
		os << endl;
	}
	
}


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
				 int NK = atoi(record.c[4]);
				for (int i = 0; i < NK; i++)
				{
					auto product = new mcEndfProduct();
					product->Load(isEndf);
					Products.push_back(product);
				}
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
	os << "Dump EA proton crossections for element = \t" << ElementName << " with \t" << Products.size() << " products." << endl;
	os << "---------------------------------------------------------------" << endl;
	os << endl;
	for (int i = 0; i < Products.size(); i++)
	{
		os << "Product - \t" << Products[i]->typeof(Products[i]->product_type) << "\t #" << i + 1 << endl;
		if (Products[i]->product_type == 5)
			os << "Nucleous with:" << endl << "A = \t" << Products[i]->ZAP % 1000 << endl << "Z = \t" << Products[i]->ZAP / 1000 << endl << endl;
		Products[i]->EANuclearCrossSections[0]->dump(os);	
	}
}


