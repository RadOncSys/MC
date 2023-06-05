#include "mcEndfP.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <random>

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

void mcEndfEANuclearCrossSectionTable::Load(std::istream& is, int LAW)
{
	string line;
	mcEndfRecord record;
	int pointCount = 0, i = 0;
	bool is_ea_read = false;
	getline(is, line, '\n');
	::memcpy(&record, line.c_str(), 80);
	if (LAW == 1)
		LANG = atoi(record.c[2]);

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
		else throw exception(("This LANG type is abcent for LAW = " + std::to_string(LAW) + ".").c_str());
	}
	else if (LAW == 2)
	{
		//counter
		int c = 0;
		EA_par[c].resize(3);
		while (!is.fail())
		{
			getline(is, line, '\n');
			::memcpy(&record, line.c_str(), 80);
			LANG = atoi(record.c[2]);
			if (LANG != 12)
				throw exception(("This LANG type is abcent for LAW = " + std::to_string(LAW) + ".").c_str());
			
			npoints_out = mcEndfRecord::iStrCrop(record.c[5], 11);
			EA_par[c].resize(npoints_out);
			double incEn = mcEndfRecord::ParseValue(record.c[1], 11);
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
	for (int i = 0; i < EANuclearCrossSections.size(); i++)
	{
		delete EANuclearCrossSections[i];
	}
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
					product_type = particle_type::gamma;
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

double mcEndfCrossSectionTable::get_sigma(double kE)
{
	int i = 0;
	for (i = 0; i < Energies.size(); i++)
		if (Energies[i] > kE)
			break;
	i--;
	return Values[i]; //сделать интерполяцию
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

double mcEndfEANuclearCrossSectionTable::playmu(double kE, double** pars, int ptype, mcRng& rng)
{
	int pi = 0;
	double SIGN = rng.rnd();
	double mur = rng.rnd();
	double r = rng.rnd();
	const double C1 = 0.04, C2 = 0.0000018, C3 = 0.00000067,
		Et1 = 130/* MeV*/, Et3 = 41/*MeV*/;
	double Ma = 0, mb = 0;

	
	if (LANG == 1) //LEGANDRE REPRESENTATION
	{
		if (SIGN > 0.5)
			return mur;
		else return -mur;
	}
	else if (LANG == 2) //KALLBACH-MANN REPRESENTATION
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
			Ap = 1;
			Zp = 0;
			break;
		case 1: mb = 1;
			Ap = 1;
			Zp = 1;
			break;
		case 2: mb = 1;
			Zp = 1;
			Ap = 2;
			break;
		case 3: mb = 1;
			Zp = 1;
			Ap = 3;
			break;
		case 4: mb = 2;
			Zp = 2;
			Ap = 4;
			break;
		case 5: throw exception("Recoils doesn't supported. See ENDF manual p.139");
			break;
		case 6: throw exception("Gammas is using another code. See ENDF manual p.139");
			break;
		default: throw exception("Unknown product particle. See ENDF manual p.139");
			break;
		}

		epsa = kE / 1000000 * AWR_nucl / (AWR_nucl + 0.99862); //Need to use AWR of incident particle instead of "0.99862"
		int AA = ZA_nucl % 1000, AC = ZA_nucl % 1000 + 1, AB = AC - Ap;
		int ZA = ZA_nucl / 1000, ZC = ZA / 1000 + 1, ZB = ZC - Zp;
		int NA = AA - ZA, NC = AC - ZC, NB = AB - ZB;
		Sa = 15.68 * (AC - AA) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NA - ZA) * (NA - ZA) / AA) - 18.56 * (pow(AC, 2 / 3) - pow(AA, 2 / 3)) +
			33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, 4 / 3) - (NA - ZA) * (NA - ZA) / pow(AA, 4 / 3)) - 0.717 * (ZC * ZC / pow(AC, 1 / 3) - ZA * ZA / pow(AA, 1 / 3)) +
			1.211 * (ZC * ZC / AC - ZA * ZA / AA); // - Ia	

		epsb = pars[0][0] / 1000000 * (AWR_nucl + 0.99862) / (AWR_nucl + 0.99862 - Ap); // AWRB = AWRA + AWRa - AWRb; in formula (p.138): AWRB + AWRb = AWRA + AWR(proton)
		Sb = 15.68 * (AC - AB) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NB - ZB) * (NB - ZB) / AB) - 18.56 * (pow(AC, 2 / 3) - pow(AB, 2 / 3)) +
			33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, 4 / 3) - (NB - ZB) * (NB - ZB) / pow(AB, 4 / 3)) - 0.717 * (ZC * ZC / pow(AC, 1 / 3) - ZB * ZB / pow(AB, 1 / 3)) +
			1.211 * (ZC * ZC / AC - ZB * ZB / AB); // - Ib

		ea = epsa + Sa;
		eb = epsb + Sb;

		double R1 = min(ea, Et1);
		double R3 = min(ea, Et3);
		a = C1 * R1 * eb / ea + C2 * pow(R1 * eb / ea, 3) + C3 * pow(R3 * eb / ea, 4) * Ma * mb;

		vector<double> mu;
		vector<double> f;
		f.resize(41);

		for (int i = -20; i < 21; i++)
		{
			mu.push_back(i * 0.05);
			if (i == -20)
				f[i + 20] = a * pars[1][0] / 2 / sinh(a) * (cosh(a * mu[i + 20]) + pars[2][0] * sinh(a * mu[i + 20]));
			else f[i + 20] = f[i + 19] + a * pars[1][0] / 2 / sinh(a) * (cosh(a * mu[i + 20]) + pars[2][0] * sinh(a * mu[i + 20]));
		}
		for (int i = 0; i < f.size(); i++)
		{
			f[i] /= f[f.size() - 1];
		}
		for (pi = 0; pi < f.size(); pi++)
		{
			if (f[pi] > r)
				break;
		}
		double output;
		if (pi == 0)
			output = -1;
		else if (pi >= mu.size())
			output = 1;
		else output = mu[pi - 1] + (mu[pi] - mu[pi - 1]) / (f[pi] - f[pi - 1]) * (r - f[pi - 1]);
		return output;
	}
}

double** mcEndfEANuclearCrossSectionTable::playpar(mcRng& rng, double kE)
{
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
	
	if (LANG == 2) //KALLBACH-MANN REPRESENTATION
	{
		double** output = new double*[3];
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
					if (ii == 0)
						f_0[ii][0] = EA_par[i][ii][1] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (getf_0(i + 1, Eout[ii]) - EA_par[i][ii][1]);
					else f_0[ii][0] = f_0[ii - 1][0] + EA_par[i][ii][1] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (getf_0(i + 1, EA_par[i][ii][0]) - EA_par[i][ii][1]);
				}
				else
				{
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * getf_0(i + 1, Eout[ii]);
					f_0[ii][0] = f_0[ii - 1][0] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * getf_0(i + 1, Eout[ii]);
				}
			}
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
					if (ii == 0)
						f_0[ii][0] = getf_0(i, Eout[ii]) + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (EA_par[i + 1][ii][1] - getf_0(i, Eout[ii]));
					else f_0[ii][0] = f_0[ii - 1][0] + getf_0(i, Eout[ii]) + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (EA_par[i + 1][ii][1] - getf_0(i, Eout[ii]));
				}
				else
				{
					f_0[ii].resize(4);
					f_0[ii][1] = Eout[ii];
					f_0[ii][3] = r_par[ii];
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * EA_par[i + 1][ii][1];
					f_0[ii][0] = f_0[ii - 1][0] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * EA_par[i + 1][ii][1];
				}
			}
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
	else if (LANG == 1) //LEGANDRE REPRESENTATION
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
					if (ii == 0)
						f_0[ii][0] = EA_par[i][ii][1] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (getf_0(i + 1, Eout[ii]) - EA_par[i][ii][1]);
					else f_0[ii][0] = f_0[ii - 1][0] + EA_par[i][ii][1] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (getf_0(i + 1, EA_par[i][ii][0]) - EA_par[i][ii][1]);
				}
				else
				{
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * getf_0(i + 1, Eout[ii]);
					f_0[ii][0] = f_0[ii - 1][0] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * getf_0(i + 1, Eout[ii]);
				}
			}
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
					if (ii == 0)
						f_0[ii][0] = getf_0(i, Eout[ii]) + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (EA_par[i + 1][ii][1] - getf_0(i, Eout[ii]));
					else f_0[ii][0] = f_0[ii - 1][0] + getf_0(i, Eout[ii]) + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * (EA_par[i + 1][ii][1] - getf_0(i, Eout[ii]));
				}
				else
				{
					f_0[ii].resize(3);
					f_0[ii][1] = Eout[ii];
					f_0[ii][2] = (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * EA_par[i + 1][ii][1];
					f_0[ii][0] = f_0[ii - 1][0] + (kE - Energies[i]) / (Energies[i + 1] - Energies[i]) * EA_par[i + 1][ii][1];
				}
			}
			for (int ii = 0; ii < f_0.size(); ii++)
				f_0[ii][0] /= f_0[f_0.size() - 1][0];
			for (pi = 0; pi < f_0.size(); pi++)
				if (f_0[pi][0] > r)
					break;
		}

		if (pi == 0) 
			for (int iii = 0; iii < 1; iii++)
				output[iii][0] = 0;
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

double mcEndfEANuclearCrossSectionTable::getMulti(double kE)
{
	int i = 0;
	for (i = 0; i < Energies.size(); i++)
		if (Energies[i] > kE)
			break;
	double multiplicity = 0;
	if (i = Energies.size())
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
	return multiplicity;
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

		// Сечения (p,n) MT = 50 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '0')
		{
			Neutron0CrossSection.Load(isEndf);
		}

		// Сечения (p,n) MT = 51 
		else if (record.MF[0] == ' ' && record.MF[1] == '3' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '1')
		{
			Neutron1CrossSection.Load(isEndf);
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
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '0')  // (p,n) реакции MT = 50, 51
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
					EmittedNeutrons.push_back(neutron);
					EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
					EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
				}
			}
		}

		else if (record.MF[0] == ' ' && record.MF[1] == '6' &&
			record.MT[0] == ' ' && record.MT[1] == '5' && record.MT[2] == '1')  // (p,n) реакции MT = 50, 51
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
					EmittedNeutrons.push_back(neutron);
					EmittedNeutrons[i]->EANuclearCrossSections[0]->AWR_nucl = AWR;
					EmittedNeutrons[i]->EANuclearCrossSections[0]->ZA_nucl = ZA;
				}
			}
		}

		getline(isEndf, line, '\n');
	}
	cout << "END";
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
		os << "Product - \t" << typeof(Products[i]->product_type) << "\t #" << i + 1 << endl;
		if (Products[i]->product_type == 5)
			os << "Nucleous with:" << endl << "A = \t" << Products[i]->ZAP % 1000 << endl << "Z = \t" << Products[i]->ZAP / 1000 << endl << endl;
		Products[i]->EANuclearCrossSections[0]->dump(os);	
	}
}


