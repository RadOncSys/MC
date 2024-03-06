#include "mcEndfP.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <random>
#include "mcSpline.h"
#include <cmath>

using namespace std;

double Legandre(int NL, double* a, double mu);

double mcEndfRecord::ParseValue(const char* s, int n)
{
	std::string s1 = s;
	s1.erase(n);
	double f = stof(s1);
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

void mcEndfEANuclearCrossSectionTable::Load(std::istream& is, int LAW)
{
	string line;
	mcEndfRecord record;
	int pointCount = 0, i = 0;
	bool is_ea_read = false;
	getline(is, line, '\n');
	::memcpy(&record, line.c_str(), 80);
	if (LAW == 1)
		LANG.push_back(atoi(record.c[2]));

	//чтение энерго-угловых параметров
	getline(is, line, '\n');
	::memcpy(&record, line.c_str(), 80);
	if (LAW == 2)
	{
		n_energypoints = atoi(record.c[0]);
		EA_par.resize(n_energypoints);
	}
	if (LAW == 1)
	{
		if (LANG[0] == 2)
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
					EA_Epoints.resize(n_energypoints);
					getline(is, line, '\n');
					is_ea_read = true;
				}
				else
				{
					for (i = 0; i < n_energypoints; i++)
					{
						NA = atoi(record.c[3]);
						EA_Epoints[i] = mcEndfRecord::ParseValue(record.c[1], 11);
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
		else if (LANG[0] == 1)
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
					EA_Epoints.resize(n_energypoints);
					getline(is, line, '\n');
					is_ea_read = true;
				}
				else
				{
					for (i = 0; i < n_energypoints; i++)
					{
						int c = 0, c1 = 0; //counters
						NA = atoi(record.c[3]);
						EA_Epoints[i] = mcEndfRecord::ParseValue(record.c[1], 11);
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
								if (c % 5 == 0 && j != npoints_out - 1)
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
		else throw exception(("This LANG = " + std::to_string(LANG[0]) + " type is abcent for LAW = " + std::to_string(LAW) + ".").c_str());
	}
	else if (LAW == 2)
	{
		//counter
		int c = 0;
		while (!is.fail())
		{
			getline(is, line, '\n');
			::memcpy(&record, line.c_str(), 80);
			LANG.push_back(atoi(record.c[2]));
			if (LANG[c] == 14)
				throw exception(("This LANG = " + std::to_string(LANG[c]) + " type is abcent for LAW = " + std::to_string(LAW) + ".").c_str());
			
			npoints_out = mcEndfRecord::iStrCrop(record.c[5], 11);
			if (LANG[c] == 12)
				EA_par[c].resize(npoints_out);
			else if (LANG[c] == 0)
				EA_par[c].resize(npoints_out + 1);
			double incEn = mcEndfRecord::ParseValue(record.c[1], 11);
			if (LANG[c] == 12)
			{
				for (int jj = 0; jj < npoints_out; jj++)
				{
					if (jj % 3 == 0)
					{
						getline(is, line, '\n');
						::memcpy(&record, line.c_str(), 80);
					}
					EA_par[c][jj].push_back(incEn);
					EA_par[c][jj].push_back(mcEndfRecord::ParseValue(record.c[(2 * jj) % 6], 11));
					EA_par[c][jj].push_back(mcEndfRecord::ParseValue(record.c[(2 * jj + 1) % 6], 11));
				}
				c++;
				if (c == n_energypoints)
					break;
			}
			else if (LANG[c] == 0)
			{
				EA_par[c][0].push_back(incEn);
				for (int jj = 0; jj < EA_par[c].size() - 1; jj++)
				{
					if (jj % 6 == 0)
					{
						getline(is, line, '\n');
						::memcpy(&record, line.c_str(), 80);
					}
					EA_par[c][jj + 1].push_back(mcEndfRecord::ParseValue(record.c[jj % 6], 11));
				}
				c++;
				if (c == n_energypoints)
					break;
			}
		}
		iLang++;
	}
	else if (LAW == 4)
	{
		getline(is, line, '\n');
		::memcpy(&record, line.c_str(), 80);
	}
	else throw exception("This LAW type doesn't supported by this code");
}

mcEndfProduct::mcEndfProduct()
{
}

mcEndfProduct::~mcEndfProduct()
{
	/*for (int i = 0; i < EANuclearCrossSections.size(); i++)
	{
		delete EANuclearCrossSections[i];
	}*/
}

std::string typeof(int i)
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
			ZAP = round(mcEndfRecord::ParseValue(record.c[0], 11));
			AWP = mcEndfRecord::ParseValue(record.c[1], 11);
			LAW = mcEndfRecord::iStrCrop(record.c[3], 11);
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
					product_type = particle_type::gammas;
					break;
				case 1000:    //dont know real electron ZAP
				case -1000:
					product_type = particle_type::electron;
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
			if (LAW != 4)
				energyangle->Load(is, LAW);
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

double mcEndfCrossSectionTable::get_lambda(double kE, double rho, double A)
{
	const double Na = 6.022; //multiply by 10^23 was taken into account, when convert barns into cm^2 at "return"
	int i = 0;
	for (i = 0; i < Energies.size(); i++)
		if (Energies[i] > kE)
			break;
	if (i == Energies.size())
		i--;
	i--;
	double _sigma = (Values[i + 1] - Values[i]) / (Energies[i + 1] - Energies[i]) * (kE - Energies[i]) + Values[i];
	return 1 / (rho * _sigma * Na / A / 10);
}

double mcEndfCrossSectionTable::get_sigma(double kE) const
{
	int i = 0;
	for (i = 0; i < Energies.size(); i++)
		if (Energies[i] > kE)
			break;
	if (i == Energies.size())
		i--;
	i--;
	double _sigma = (Values[i + 1] - Values[i]) / (Energies[i + 1] - Energies[i]) * (kE - Energies[i]) + Values[i];
	return _sigma;
}

void mcEndfEANuclearCrossSectionTable::dump(std::ostream& os) const
{
	os << "Multiplicity table:" << endl;
	os << "# \t" << "Energy \t" << "Multiplicity" << endl;
	for (int i = 0; i < n_energypoints; i++)
	{
		os << i << "\t" << Energies[i] << "\t" << Multiplicities[i] << endl;
	}
	if (LANG[0] == 2)
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
	else if (LANG[0] == 1)
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

double mcEndfEANuclearCrossSectionTable::playmu(double kE, int LAW, int keIN, int eoutID, int ptype, mcRng& rng) const
{
	int pi = 0;
	double SIGN = rng.rnd();
	double mur = rng.rnd();
	double r = rng.rnd();
	double ksi1 = rng.rnd();
	double ksi2 = rng.rnd();
	const double C1 = 0.04, C2 = 0.0000018, C3 = 0.00000067,
		Et1 = 130/* MeV*/, Et3 = 41/*MeV*/;
	double Ma = 0, mb = 0;
	if (LAW == 1)
	{
		if (LANG[0] == 1) //LEGANDRE REPRESENTATION
		{
			if (SIGN > 0.5)
				return mur;
			else return -mur;
		}
		else if (LANG[0] == 2) //KALLBACH-MANN REPRESENTATION
		{
			double a = 0; // kallbach-mann parameter
			double epsa = 0, epsb = 0, Sa = 0, Sb = 0;
			double ea = 0, eb = 0;
			int Zp, Ap;

			switch (1) {								//Need to use type of incedent particle instead of "1"
			case 0: Ma = 1;
				break;
			case 1: Ma = 1;
				break;
			case 2: Ma = 1;
				break;
			case 3: throw exception("Triton doesn't supported. See ENDF manual p.139");
				break;
			case 4: Ma = 0;
				break;
			case 5: throw exception("Recoils doesn't supported. See ENDF manual p.139");
				break;
			case 6: throw exception("Gammas is using another code. See ENDF manual p.139");
				break;
			default: throw exception("Unknown product particle. See ENDF manual p.139");
				break;
			}
			switch (ptype) {
			case 0: mb = 0.5;
				Ap = 1.0;
				Zp = 0.0;
				break;
			case 1: mb = 1.0;
				Ap = 1.0;
				Zp = 1.0;
				break;
			case 2: mb = 1.0;
				Zp = 1.0;
				Ap = 2.0;
				break;
			case 3: mb = 1.0;
				Zp = 1.0;
				Ap = 3.0;
				break;
			case 4: mb = 2.0;
				Zp = 2.0;
				Ap = 4.0;
				break;
			case 5: throw exception("Recoils doesn't supported. See ENDF manual p.139");
				break;
			case 6: throw exception("Gammas is using another code. See ENDF manual p.139");
				break;
			default: throw exception("Unknown product particle. See ENDF manual p.139");
				break;
			}

			epsa = kE / 1000000 * AWR_nucl / (AWR_nucl + 0.99862); //Need to use AWR of incident particle instead of "0.99862"
			double AA = ZA_nucl % 1000, AC = ZA_nucl % 1000 + 1, AB = AC - Ap;
			double ZA = ZA_nucl / 1000, ZC = ZA + 1, ZB = ZC - Zp;
			double NA = AA - ZA, NC = AC - ZC, NB = AB - ZB;
			double kf1 = 4.0 / 3.0, kf2 = 2.0 / 3.0, kf3 = 1.0 / 3.0;
			Sa = 15.68 * (AC - AA) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NA - ZA) * (NA - ZA) / AA) - 18.56 * (pow(AC, kf2) - pow(AA, kf2)) +
				33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, kf1) - (NA - ZA) * (NA - ZA) / pow(AA, kf1)) - 0.717 * (ZC * ZC / pow(AC, kf3) - ZA * ZA / pow(AA, kf3)) +
				1.211 * (ZC * ZC / AC - ZA * ZA / AA); // - Ia	

			epsb = EA_par[keIN][eoutID][0] / 1000000 * (AWR_nucl + 0.99862) / (AWR_nucl + 0.99862 - Ap); // AWRB = AWRA + AWRa - AWRb; in formula (p.138): AWRB + AWRb = AWRA + AWR(proton)
			Sb = 15.68 * (AC - AB) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NB - ZB) * (NB - ZB) / AB) - 18.56 * (pow(AC, kf2) - pow(AB, kf2)) +
				33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, kf1) - (NB - ZB) * (NB - ZB) / pow(AB, kf1)) - 0.717 * (ZC * ZC / pow(AC, kf3) - ZB * ZB / pow(AB, kf3)) +
				1.211 * (ZC * ZC / AC - ZB * ZB / AB); // - Ib

			ea = epsa + Sa;
			eb = epsb + Sb;

			double R1 = min(ea, Et1);
			double R3 = min(ea, Et3);
			a = C1 * R1 * eb / ea + C2 * pow(R1 * eb / ea, 3) + C3 * pow(R3 * eb / ea, 4) * Ma * mb;

			double Excl = (1 + EA_par[keIN][eoutID][2]) / 2;
			double mu = 0;
			if (ksi1 < Excl)
				mu = 1 / a * log(1 + ksi2 * (exp(2 * a) - 1)) - 1;
			else
				mu = -1 / a * log(1 + ksi2 * (exp(-2 * a) - 1)) - 1;
			return mu;
		}
	}
	else if (LAW == 2 && EA_par[keIN][eoutID][0])
	{
		throw exception ("Trying to play particle direction using LAW = 2.");
		int interpol = 0;
		for (interpol = 0; interpol < EA_par.size(); interpol++)
		{
			if (EA_par[interpol][0][0] > kE)
				break;
		}
		if (interpol != EA_par.size())
		{
			if (LANG[interpol] == 0)
			{
				vector<double> mu;
				vector<double> p;
				double* a;
				int SIZE = EA_par[interpol].size();
				a = new double[SIZE];
				for (int i = 0; i < EA_par[interpol].size(); i++)
				{
					a[i] = EA_par[interpol][i][0];
				}
				for (int i = -20; i < 21; i++)
					mu.push_back(i * 0.05);
				p.resize(mu.size());
				for (int i = 0; i < mu.size(); i++)
				{
					if (i == 0)
						p[i] = (Legandre(EA_par[interpol].size(), a, mu[i]));
					else p[i] = p[i - 1] + (Legandre(EA_par[interpol].size(), a, mu[i]));
				}
				for (int i = 0; i < mu.size(); i++)
				{
					p[i] /= p[p.size() - 1];
				}
			}
		}
		else 
		{

		}
	}
}

double Pn(double x, int n)
{
	if (n == 0) return 1;
	if (n == 1) return x;
	return ((2 * n - 1) * Pn(x, n - 1) - (n - 1) * Pn(x, n - 2)) / n;
}

double Legandre(int NL, double* a, double mu)
{
	double output = 0.5;
	for (double j = 1.0; j < NL; j++)
		output += legendre(int(j), mu) * a[int(j)] * (2.0 * j + 1.0) / 2.0;
	return output;
}

double LegandreKM(int NA, double* f, double mu)
{
	double output = 0;
	for (double j = 1.0; j < NA; j++)
		output += pow(mu, j) * f[int(j)] * (2.0 * j + 1.0) / 2.0;
	return output;
}


double** mcEndfEANuclearCrossSectionTable::playpar(mcRng& rng, double kE, int LAW)
{
	if (LAW == 2)
	{
		double** output = new double* [1];
		for (int i = 0; i < 1; i++)
		{
			output[i] = new double[1];
		}
		output[0][0] = 0;
		return output;
	}
	//f(e, e') = f(ei, e') + (e - ei)/(ei+1 - ei)*(f(ei+1, e') - f(ei,e'))
	int i = 0, c = 0, iii = 0;
	int pi = 0;
	//double r = double(rand()) / RAND_MAX; 
	vector<vector<double>> f_0;
	double r = rng.rnd();
	for (i = 0; i < EA_Epoints.size(); i++)
	{
		if (EA_Epoints[i] > kE)
			break;
	}
	if (i == EA_Epoints.size())
		i -= 2;
	else i--;
	int maxsize = max(EA_par[i].size(), EA_par[i + 1].size());
	int minsize = min(EA_par[i].size(), EA_par[i + 1].size());
	vector<double> Eout;
	vector<double> r_par;
	//vector<vector<double>> f_l (2);
	
	if (LANG[0] == 2) //KALLBACH-MANN REPRESENTATION
	{
		double** output = new double* [3];
		for (int i = 0; i < 3; i++)
		{
			output[i] = new double[1];
		}

		if (EA_par[i].size() > EA_par[i + 1].size())
		{
			for (int ii = 0; ii < maxsize; ii++)
			{
				Eout.push_back(EA_par[i][ii][0]);
				r_par.push_back(EA_par[i][ii][2]);
			}
			for (c = 0; c < minsize; c++)
			{
				if (EA_par[i + 1][c][0] > Eout.back())
					break;
			}
			for (iii = c; iii < minsize; iii++)
			{
				Eout.push_back(EA_par[i + 1][iii][0]);
				r_par.push_back(EA_par[i + 1][iii][2]);
			}
			f_0.resize(Eout.size());
			for (int ii = 0; ii < Eout.size(); ii++)
			{
				f_0[ii].resize(4);
				f_0[ii][1] = Eout[ii];
				f_0[ii][3] = r_par[ii];
				if (ii < maxsize)
				{
					f_0[ii][2] = EA_par[i][ii][1] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (getf_0(i + 1, Eout[ii]) - EA_par[i][ii][1]);
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
				else
				{
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * getf_0(i + 1, Eout[ii]);
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
			}
			f_0[f_0.size() - 1][0] = 0;
			for (int ii = 1; ii < f_0.size(); ii++)
				f_0[ii][0] += f_0[ii - 1][0];
			for (int ii = 0; ii < f_0.size(); ii++)
				f_0[ii][0] /= f_0[f_0.size() - 1][0];
			for (pi = 0; pi < f_0.size(); pi++)
				if (f_0[pi][0] > r)
					break;
		}
		else if (EA_par[i].size() <= EA_par[i + 1].size())						
		{																		
			for (int ii = 0; ii < maxsize; ii++)
			{
				Eout.push_back(EA_par[i + 1][ii][0]);
				r_par.push_back(EA_par[i + 1][ii][2]);
			}
			f_0.resize(Eout.size());
			for (int ii = 0; ii < Eout.size(); ii++)
			{
				if (ii < minsize)
				{
					f_0[ii].resize(4);
					f_0[ii][1] = Eout[ii];
					f_0[ii][2] = getf_0(i, Eout[ii]) + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (EA_par[i + 1][ii][1] - getf_0(i, Eout[ii]));
					f_0[ii][3] = r_par[ii];
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
				else
				{
					f_0[ii].resize(4);
					f_0[ii][1] = Eout[ii];
					f_0[ii][3] = r_par[ii];
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * EA_par[i + 1][ii][1];
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
			}
			f_0[f_0.size() - 1][0] = 0;
			for (int ii = 1; ii < f_0.size(); ii++)
				f_0[ii][0] += f_0[ii - 1][0];
			for (int ii = 0; ii < f_0.size(); ii++)
				f_0[ii][0] /= f_0[f_0.size() - 1][0];
			for (pi = 0; pi < f_0.size(); pi++)
				if (f_0[pi][0] > r)
					break;
		}
		if (pi == 0) for (int iii = 0; iii < 3; iii++)
			output[iii][0] = 0;
		else if (pi >= f_0.size())
		{
			output[0][0] = f_0[f_0.size() - 1][1];
			output[1][0] = f_0[f_0.size() - 1][2];
			output[2][0] = f_0[f_0.size() - 1][3];
		}
		else if (pi != 0)
		{
			output[0][0] = f_0[pi - 1][1] + (f_0[pi][1] - f_0[pi - 1][1]) / (f_0[pi][0] - f_0[pi - 1][0]) * (r - f_0[pi - 1][0]);	//Eout
			output[1][0] = f_0[pi - 1][2] + (f_0[pi][2] - f_0[pi - 1][2]) / (f_0[pi][0] - f_0[pi - 1][0]) * (r - f_0[pi - 1][0]);	//f_0
			output[2][0] = f_0[pi - 1][3] + (f_0[pi][3] - f_0[pi - 1][3]) / (f_0[pi][0] - f_0[pi - 1][0]) * (r - f_0[pi - 1][0]);	//r
		}
		return output;
	}
	else if (LANG[0] == 1) //LEGANDRE REPRESENTATION
	{
		bool imax = true;
		double** output = new double* [2];
		output[0] = new double[1];
		output[1] = new double[1]; //if NA > 0 expand this array and use vector f_l (declared in the beginning of the fuction) to make Legandre polinom

		if (EA_par[i].size() > EA_par[i + 1].size())
		{
			for (int ii = 0; ii < maxsize; ii++)
			{
				Eout.push_back(EA_par[i][ii][0]);
			}
			for (c = 0; c < minsize; c++)
			{
				if (EA_par[i + 1][c][0] > Eout.back())
					break;
			}
			for (iii = c; iii < minsize; iii++)
			{
				Eout.push_back(EA_par[i + 1][iii][0]);
			}
			f_0.resize(Eout.size());
			for (int ii = 0; ii < Eout.size(); ii++)
			{
				f_0[ii].resize(3);
				f_0[ii][1] = Eout[ii];
				if (ii < maxsize)
				{
					f_0[ii][2] = EA_par[i][ii][1] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (getf_0(i + 1, Eout[ii]) - EA_par[i][ii][1]);
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
				else
				{
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * getf_0(i + 1, Eout[ii]);
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
			}
			f_0[f_0.size() - 1][0] = 0;
			for (int ii = 1; ii < f_0.size(); ii++)
				f_0[ii][0] += f_0[ii - 1][0];
			for (int ii = 0; ii < f_0.size(); ii++)
				f_0[ii][0] /= f_0[f_0.size() - 1][0];
			for (pi = 0; pi < f_0.size(); pi++)
				if (f_0[pi][0] > r)
					break;
		}
		else if (EA_par[i].size() <= EA_par[i + 1].size())
		{
			imax = false;
			for (int ii = 0; ii < maxsize; ii++)
			{
				Eout.push_back(EA_par[i + 1][ii][0]);
			}
			f_0.resize(Eout.size());
			for (int ii = 0; ii < Eout.size(); ii++)
			{
				if (ii < minsize)
				{
					f_0[ii].resize(3);
					f_0[ii][1] = Eout[ii];
					f_0[ii][2] = getf_0(i, Eout[ii]) + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (EA_par[i + 1][ii][1] - getf_0(i, Eout[ii]));
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
				else
				{
					f_0[ii].resize(3);
					f_0[ii][1] = Eout[ii];
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * EA_par[i + 1][ii][1];
					if (ii < Eout.size() - 1)
						f_0[ii][0] = f_0[ii][2] * (Eout[ii + 1] - Eout[ii]);
				}
			}
			f_0[f_0.size() - 1][0] = 0;
			for (int ii = 1; ii < f_0.size(); ii++)
				f_0[ii][0] += f_0[ii - 1][0];
			for (int ii = 0; ii < f_0.size(); ii++)
				f_0[ii][0] /= f_0[f_0.size() - 1][0];
			for (pi = 0; pi < f_0.size(); pi++)
				if (f_0[pi][0] > r)
					break;
		}

		if (pi == 0)
		{
			output[0][0] = 0;
			output[1][0] = f_0[0][0];
		}
		else if (pi >= f_0.size())
		{
			output[0][0] = f_0[f_0.size() - 1][1];
			output[1][0] = f_0[f_0.size() - 1][2];
		}
		else if (pi != 0)
		{
			output[0][0] = f_0[pi - 1][1] + (f_0[pi][1] - f_0[pi - 1][1]) / (f_0[pi][0] - f_0[pi - 1][0]) * (r - f_0[pi - 1][0]);
			output[1][0] = f_0[pi - 1][2] + (f_0[pi][2] - f_0[pi - 1][2]) / (f_0[pi][0] - f_0[pi - 1][0]) * (r - f_0[pi - 1][0]);
		}
		return output;
	}
}

double lagrange(double xp, vector<vector<double>> f) {
	double intp = 0;
	int n = f.size();
	for (int i = 0; i < n; i++) {
		double m = 1;
		for (int j = 0; j < n; j++) {
			if (i != j)
				m = m * (xp - f[j][1]) / (f[i][1] - f[j][1]);
		}
		m = m * f[i][0];
		intp = intp + m;
	}
	return intp;
}

double LinearInt(double X, vector<vector<double>> f)
{
	int i = 0;
	for (i = 0; i < f.size(); i++)
		if (f[i][1] > X)
			break;
	if (i >= f.size())
		i -= 2;
	else i--;
	double output = f[i][0] + (f[i + 1][0] - f[i][0]) / (f[i + 1][1] - f[i][1]) * (X - f[i][1]);
	return output;
}

double mcEndfEANuclearCrossSectionTable::integrate_f0(mcRng& rng, double kE)
{
	//f(e, e') = f(ei, e') + (e - ei)/(ei+1 - ei)*(f(ei+1, e') - f(ei,e'))
	int i = 0, c = 0, iii = 0;
	int pi = 0;
	vector<vector<double>> f_0;
	double r = rng.rnd();
	for (i = 0; i < EA_Epoints.size(); i++)
	{
		if (EA_Epoints[i] > kE)
			break;
	}
	if (i == EA_Epoints.size())
		i -= 2;
	else i--;
	int maxsize = max(EA_par[i].size(), EA_par[i + 1].size());
	int minsize = min(EA_par[i].size(), EA_par[i + 1].size());
	vector<double> Eout;
	vector<double> r_par;
	//vector<vector<double>> f_l (2);

	if (LANG[0] == 2) //KALLBACH-MANN REPRESENTATION
	{
		if (EA_par[i].size() > EA_par[i + 1].size())
		{
			for (int ii = 0; ii < maxsize; ii++)
			{
				Eout.push_back(EA_par[i][ii][0]);
			}
			for (c = 0; c < minsize; c++)
			{
				if (EA_par[i + 1][c][0] > Eout.back())
					break;
			}
			for (iii = c; iii < minsize; iii++)
			{
				Eout.push_back(EA_par[i + 1][iii][0]);
			}
			f_0.resize(Eout.size());
			for (int ii = 0; ii < Eout.size(); ii++)
			{
				f_0[ii].resize(2);
				f_0[ii][1] = Eout[ii];
				if (ii < maxsize)
				{
					f_0[ii][0] = EA_par[i][ii][1] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (getf_0(i + 1, Eout[ii]) - EA_par[i][ii][1]);
				}
				else
				{
					f_0[ii][0] = f_0[ii - 1][0] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * getf_0(i + 1, Eout[ii]);
				}
			}
		}
		else if (EA_par[i].size() <= EA_par[i + 1].size())
		{
			for (int ii = 0; ii < maxsize; ii++)
			{
				Eout.push_back(EA_par[i + 1][ii][0]);
			}
			f_0.resize(Eout.size());
			for (int ii = 0; ii < Eout.size(); ii++)
			{
				if (ii < minsize)
				{
					f_0[ii].resize(2);
					f_0[ii][1] = Eout[ii];
					f_0[ii][0] = getf_0(i, Eout[ii]) + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (EA_par[i + 1][ii][1] - getf_0(i, Eout[ii]));
				}
				else
				{
					f_0[ii].resize(2);
					f_0[ii][1] = Eout[ii];
					f_0[ii][0] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * EA_par[i + 1][ii][1];
				}
			}
		}
		
	}
	vector<double> E, F;
	for (int i = 0; i < f_0.size(); i++)
	{
		E.push_back(f_0[i][1]);
		F.push_back(f_0[i][0]);
	}
	vector<SplineSet> cs = spline(E, F);
	double S = 0;
	vector<vector<double>> f1;
	f1.resize(1000);
	double Eh = 0;
	double h = f_0[f_0.size() - 1][1] / f1.size();
	for (int i = 0; i < f1.size(); i++)
	{
		f1[i].resize(2);
		f1[i][1] = Eh;
		f1[i][0] = CountSpline(cs, f1[i][1]);
		Eh += h;
		//cout << f1[i][1] << "\t" << f1[i][0] << endl;
	}
	for (int i = 1; i < f1.size() - 1; i += 2)
	{
		S += h / 3 * (f1[i - 1][0] + 4 * f1[i][0] + f1[i + 1][0]);
	}
	return S;
}

double mcEndfEANuclearCrossSectionTable::getf_0(int IN, double Eout)
{
	int i = 0;
	for (i = 0; i < EA_par[IN].size(); i++)
		if (EA_par[IN][i][0] > Eout)
			break;
	if (i == EA_par[IN].size())
		i--;
	double f_0 = 0;
	if (i > 0)
		switch (interpolationType)
		{
		case 2: //lin-lin
			f_0 = (Eout - EA_par[IN][i - 1][0]) / (EA_par[IN][i][0] - EA_par[IN][i - 1][0]) * (EA_par[IN][i][1] - EA_par[IN][i - 1][1]) + EA_par[IN][i - 1][1];
			break;
		case 3: //lin-log
			f_0 = (log(Eout) - log(EA_par[IN][i][0])) / (log(EA_par[IN][i][0]) - log(EA_par[IN][i - 1][0])) * (EA_par[IN][i][1] - EA_par[IN][i - 1][1]) + EA_par[IN][i - 1][1];
			break;
		case 4: //log-log
			f_0 = exp((Eout - EA_par[IN][i][0]) / (EA_par[IN][i][0] - EA_par[IN][i - 1][0]) * (log(EA_par[IN][i][1]) - log(EA_par[IN][i - 1][1])) + log(EA_par[IN][i - 1][1]));
			break;
		default:  //undefined -> lin-lin
			f_0 = (Eout - EA_par[IN][i - 1][0]) / (EA_par[IN][i][0] - EA_par[IN][i - 1][0]) * (EA_par[IN][i][1] - EA_par[IN][i - 1][1]) + EA_par[IN][i - 1][1];
		}
	else throw exception("Interpolation range break");
	return f_0;
}

int mcEndfEANuclearCrossSectionTable::playMulti(double kE, mcRng& rng) const
{
	int i = 0;
	for (i = 0; i < Energies.size(); i++)
		if (Energies[i] > kE)
			break;
	double multiplicity = 0;
	if (i == Energies.size())
		i--;
	if (i > 0)
		switch (interpolationType)
		{
		case 2: //lin-lin
			multiplicity = (kE - Energies[i - 1]) / (Energies[i] - Energies[i - 1]) * (Multiplicities[i] - Multiplicities[i - 1]) + Multiplicities[i - 1];
			break;
		case 3: //lin-log
			multiplicity = (log(kE) - log(Energies[i - 1])) / (log(Energies[i]) - log(Energies[i - 1])) * (Multiplicities[i] - Multiplicities[i - 1]) + Multiplicities[i - 1];
			break;
		case 4: //log-log
			multiplicity = exp((kE - Energies[i - 1]) / (Energies[i] - Energies[i - 1]) * (log(Multiplicities[i]) - log(Multiplicities[i - 1])) + log(Multiplicities[i - 1]));
			break;
		default: //undefined -> lin-lin
			multiplicity = (kE - Energies[i - 1]) / (Energies[i] - Energies[i - 1]) * (Multiplicities[i] - Multiplicities[i - 1]) + Multiplicities[i - 1];
		}
	else throw exception("Interpolation range break");

	int quantity = int(multiplicity);
	double additional = multiplicity - quantity;
	double random = rng.rnd();
	if (random > additional)
		quantity++;

	return quantity;
}

double mcEndfEANuclearCrossSectionTable::playE(double kE, int &keIN, int &eoutID, mcRng& rng) const
{
	eoutID = 0;
	vector<double> probability;
	kE *= 1000000;
	keIN = 0;
	for (keIN = 0; keIN < Energies.size(); keIN++)
		if (Energies[keIN] > kE)
			break;
	if (keIN == 0)
	{
		for (int i = 0; i < EA_par[keIN].size(); i++)
		{
			if (i != EA_par[keIN].size() - 1)
				probability.push_back(EA_par[keIN][i][1] * (EA_par[keIN][i + 1][0] - EA_par[keIN][i][0]));
			else
				probability.push_back(0);
		}
	}
	else if (Energies[keIN] - kE > kE - Energies[keIN - 1])
	{
		keIN--;
		for (int i = 0; i < EA_par[keIN].size(); i++)
		{
			if (i != EA_par[keIN].size() - 1)
				probability.push_back(EA_par[keIN][i][1] * (EA_par[keIN][i + 1][0] - EA_par[keIN][i][0]));
			else
				probability.push_back(0);
		}
	}
	else
	{
		for (int i = 0; i < EA_par[keIN].size(); i++)
		{
			if (i != EA_par[keIN].size() - 1)
				probability.push_back(EA_par[keIN][i][1] * (EA_par[keIN][i + 1][0] - EA_par[keIN][i][0]));
			else
				probability.push_back(0);
		}
	}
	double random = rng.rnd();
	double psum = 0;
	for (int i = 0; i < probability.size(); i++)
		psum += probability[i];
	for (int i = 0; i < probability.size(); i++)
		probability[i] /= psum;
	for (int i = 1; i < probability.size(); i++)
	{
		probability[i] += probability[i - 1];
	}
	for (eoutID = 0; eoutID < probability.size(); eoutID++)
		if (probability[eoutID] > random)
			break;

	return EA_par[keIN][eoutID][0] / 1000000.0;
}


mcEndfP::mcEndfP()
{
}

mcEndfP::~mcEndfP()
{
	/*for (int i = 0; i < Products.size(); i++)
	{
		delete Products[i];
	}*/
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
			TotalCrossSections.isEmpty = true;
			NuclearCrossSections.isEmpty = true;
			Neutron0CrossSection.isEmpty = true;
			Neutron1CrossSection.isEmpty = true;
			Neutron2CrossSection.isEmpty = true;
			Neutron3CrossSection.isEmpty = true;
			Neutron4CrossSection.isEmpty = true;
			Neutron5CrossSection.isEmpty = true;
		
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
			TotalCrossSections.isEmpty = false;
		}

		// Сечения ядерных реакций
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')
		{
			NuclearCrossSections.Load(isEndf);
			NuclearCrossSections.isEmpty = false;
		}

		// Сечения (p,n) MT = 50 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '0')
		{
			Neutron0CrossSection.Load(isEndf);
			Neutron0CrossSection.isEmpty = false;
		}

		// Сечения (p,n) MT = 51 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '1')
		{
			Neutron1CrossSection.Load(isEndf);
			Neutron1CrossSection.isEmpty = false;
		}

		// Сечения (p,n) MT = 52 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '2')
		{
			Neutron2CrossSection.Load(isEndf);
			Neutron2CrossSection.isEmpty = false;
		}

		// Сечения (p,n) MT = 53 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '3')
		{
			Neutron3CrossSection.Load(isEndf);
			Neutron3CrossSection.isEmpty = false;
		}

		// Сечения (p,n) MT = 54 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '4')
		{
			Neutron4CrossSection.Load(isEndf);
			Neutron4CrossSection.isEmpty = false;
		}

		// Сечения (p,n) MT = 55 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '5')
		{
			Neutron5CrossSection.Load(isEndf);
			Neutron5CrossSection.isEmpty = false;
		}

		// Энерго-угловые распределения
		else if (record.MF[0] == ' ' && record.MF[1] == '6' && 
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')  //ядерные реакции (остаточные)
		{
			if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
			{
				 int NK = atoi(record.c[4]);
				 int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
				 double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
				for (int i = 0; i < NK; i++)
				{
					auto product = new mcEndfProduct();
					product->Load(isEndf);
					Products.push_back(product);
					Products[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
					Products[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
				}
			}
		}

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '0')  // (p,n) реакции MT = 50 MF = 6
		{
			if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
			{
				int NK = atoi(record.c[4]);
				int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
				double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
				for (int i = 0; i < NK; i++)
				{
					auto neutron = new mcEndfProduct();
					neutron->Load(isEndf);
					if (neutron->product_type == 0)
					{
						EmittedNeutrons.push_back(neutron);
						EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
						EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
					}
				}
			}
		}

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '1')  // (p,n) реакции MT = 51 MF = 6
		{
			if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
			{
				int NK = atoi(record.c[4]);
				int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
				double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
				for (int i = 0; i < NK; i++)
				{
					auto neutron = new mcEndfProduct();
					neutron->Load(isEndf);
					if (neutron->product_type == 0)
					{
						EmittedNeutrons.push_back(neutron);
						EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
						EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
					}
				}
			 }
		}

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '2')  // (p,n) реакции MT = 52 MF = 6
			{
				if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
				{
					int NK = atoi(record.c[4]);
					int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
					double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
					for (int i = 0; i < NK; i++)
					{
						auto neutron = new mcEndfProduct();
						neutron->Load(isEndf);
						if (neutron->product_type == 0)
						{
							EmittedNeutrons.push_back(neutron);
							EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
							EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
						}
					}
				}
			 }

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '3')  // (p,n) реакции MT = 53 MF = 6
			{
				if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
				{
					int NK = atoi(record.c[4]);
					int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
					double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
					for (int i = 0; i < NK; i++)
					{
						auto neutron = new mcEndfProduct();
						neutron->Load(isEndf);
						if (neutron->product_type == 0)
						{
							EmittedNeutrons.push_back(neutron);
							EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
							EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
						}
					}
				}
			 }

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '4')  // (p,n) реакции MT = 54 MF = 6
			{
				if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
				{
					int NK = atoi(record.c[4]);
					int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
					double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
					for (int i = 0; i < NK; i++)
					{
						auto neutron = new mcEndfProduct();
						neutron->Load(isEndf);
						if (neutron->product_type == 0)
						{
							EmittedNeutrons.push_back(neutron);
							EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
							EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
						}
					}
				}
			 }

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '5')  // (p,n) реакции MT = 55 MF = 6
			{
				if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
				{
					int NK = atoi(record.c[4]);
					int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
					double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
					for (int i = 0; i < NK; i++)
					{
						auto neutron = new mcEndfProduct();
						neutron->Load(isEndf);
						if (neutron->product_type == 0)
						{
							EmittedNeutrons.push_back(neutron);
							EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
							EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
						}
					}
				}
			 }

		getline(isEndf, line, '\n');
	}
}

void mcEndfAngular::Load(std::istream& is)
{
	string line, s1, s2, s3, s4;
	mcEndfRecord record;
	int pointCount = 0;
	bool isFirstTime = true;
	bool isNE1read = false;
	isEmpty = false;

	getline(is, line, '\n');

	while (!is.fail())
	{
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);
		if (LTT == 3)
		{
			if (isFirstTime)
			{
				LI = mcEndfRecord::ParseValue(record.c[3], 11);
				LCT = mcEndfRecord::ParseValue(record.c[4], 11);
				isFirstTime = false;
			}
			else
			{
				getline(is, line, '\n');
				::memcpy(&record, line.c_str(), 80);
				if (!isNE1read)
				{
					NE1 = mcEndfRecord::ParseValue(record.c[5], 11);
				}
			}
		}
	}
}

mcEndfN::mcEndfN()
{
}

mcEndfN::~mcEndfN()
{
}

void mcEndfN::Load(const char* fname, const char* ename)
{
	// Separators
	static const string beginSeparator("*** C O N T E N T S ***");
	static const string beginSeparator2("*****");

	ifstream isEndf(fname);
	if (isEndf.fail())
		throw exception((string("Can't open Neutron data file: ") + ename).c_str());
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
		if (line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		::memcpy(&record, line.c_str(), 80);

		// Начало новой энергии
		if (!isInData)
		{
			if (line.find(beginSeparator) != string::npos || line.find(beginSeparator2) != string::npos)
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

		// Полное сечение
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '1')
		{
			TotalCrossSections.Load(isEndf);
			TotalCrossSections.isEmpty = false;
		}

		// Сечения упругого рассеяния	
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '2')
		{
			ElasticCrossSections.Load(isEndf);
			ElasticCrossSections.isEmpty = false;
		}

		// Сечения ядерных реакций
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')
		{
			NuclearCrossSections.Load(isEndf);
			NuclearCrossSections.isEmpty = false;
		}

		// Сечения (n,n') MT = 51 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '1')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 52
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '2')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 53 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '3')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 54 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '4')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 55 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '5')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 56 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '6')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 57 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '7')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 58 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '8')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 59 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '9')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 60 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '0')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 61 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '1')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 62		
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '2')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 63 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '3')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 64 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '4')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 65
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '5')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 66
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '6')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 67
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '7')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 68
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '8')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 69 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '6' && record.MT[2] == '9')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 70 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '0')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 71 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '1')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 72 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '2')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 73 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '3')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 74 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '4')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 75 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '5')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 76 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '6')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 77 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '7')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 78 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '8')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 79 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '7' && record.MT[2] == '9')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 80 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '0')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 81 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '1')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 82 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '2')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 83 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '3')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 84 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '4')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 85 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '5')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 86 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '6')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 87 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '7')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 88 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '8')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 89 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '8' && record.MT[2] == '9')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 90 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '9' && record.MT[2] == '0')
		{
			auto inelastic = new mcEndfCrossSectionTable();
			inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
			inelastic->Load(isEndf);
			nInelasticCS.push_back(inelastic);
		}

		// Сечения (n,n') MT = 91
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '9' && record.MT[2] == '1')
			{
				auto inelastic = new mcEndfCrossSectionTable();
				inelastic->MT = (record.MT[1] - '0') * 10 + record.MT[2] - '0';
				inelastic->Load(isEndf);
				nInelasticCS.push_back(inelastic);
				}

		else if (record.MF[0] == ' ' && record.MF[1] == '4' &&
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '2')
		{
				nElasticAngular.ZA = mcEndfRecord::ParseValue(record.c[0], 11);
				nElasticAngular.AWR = mcEndfRecord::ParseValue(record.c[1], 11);
				nElasticAngular.LTT = mcEndfRecord::ParseValue(record.c[3], 11);
				nElasticAngular.Load(isEndf);
		}


		// Энерго-угловые распределения
		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')  //ядерные реакции (остаточные)
		{
			if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
			{
				int NK = atoi(record.c[4]);
				int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
				double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
				for (int i = 0; i < NK; i++)
				{
					auto product = new mcEndfProduct();
					product->Load(isEndf);
					Products.push_back(product);
					Products[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
					Products[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
				}
			}
		}

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '9' && record.MT[2] == '1')  // (n,n) реакции MT = 91 MF = 6
		{
			if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '1')
			{
				int NK = atoi(record.c[4]);
				int ZA = mcEndfRecord::ParseValue(record.c[0], 11);
				double AWR = mcEndfRecord::ParseValue(record.c[1], 11);
				for (int i = 0; i < NK; i++)
				{
					auto continuousN = new mcEndfProduct();
					continuousN->Load(isEndf);
					nInelasticContin.push_back(continuousN);
					nInelasticContin[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
					nInelasticContin[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
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

	mcRng rng1;
	rng1.init(23, 71);

	for (int i = 0; i < Products.size(); i++)
	{
		os << "Product - \t" << typeof(Products[i]->product_type) << "\t #" << i + 1 << endl;
		if (Products[i]->product_type == 5)
			os << "Nucleous with:" << endl << "A = \t" << Products[i]->ZAP % 1000 << endl << "Z = \t" << Products[i]->ZAP / 1000 << endl << endl;
		Products[i]->EANuclearCrossSections[0]->dump(os);	
	}
	double** a = Products[0]->EANuclearCrossSections[0]->playpar(rng1, 70000000, Products[0]->LAW);
}


