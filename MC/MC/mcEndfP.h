// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
// Classes to manage cross sections load for proton nuclear interractions
// from ENDF format
//---------------------------------------------------------------------------
#pragma once
#include <string>
#include <vector>
#include <memory>

// ��������� ������� ENDF �����
struct mcEndfRecord
{
	char c[6][11];
	//char data[66];
	char Z[2];
	char Stblt[2];
	char MF[2];
	char MT[3];
	char LineNumber[5];

	// ������� ������� � ��������� ������ � ������� ENDF
	// (��� ������ ������� ������������ ����� ����� +/-).
	static double ParseValue(const char* s, int n);
};

// Crossections per istope per incident particle energy
class mcEndfCrossSectionTable
{
public:
	void Load(std::istream& is);
	void dump(std::ostream& os) const;

	// ���������� ��� ������� �������� ������� / �������
	int npoints;

	// ���������� ����� ������������
	// �������� ������������, ��� �� �� ���������� �� 
	// ��������������� ������������ � �������� ����� �������.
	// TODO: ��������� ��� �������� (� ���� ����� exception �� ���� ������)
	int ninterpolations;
	 
	// ��� ������������.
	// TODO: ���� ����������� ����������� � ��������� 
	// ��������� ����� ������������ ���������� � ������
	int interpolationType;
	
	// �����
	std::vector<double> Energies;
	std::vector<double> Values;
};

// ����� ��� ������ MF=6 MT=5 

class mcEndfEANuclearCrossSectionTable
{
public:
	void mLoad(std::istream& is);
	//void dump(std::ostream& os) const;

	void Load(std::istream& is);

	// ���������� ��� ������� �������� ������� / ���������������
	int n_energypoints;

	// ���������� ����-�� ���
	int npoints_ang;

	// ���������� ����� ������������
	// �������� ������������, ��� �� �� ���������� �� 
	// ��������������� ������������ � �������� ����� �������.
	// TODO: ��������� ��� �������� (� ���� ����� exception �� ���� ������)
	int ninterpolations;

	// ��� ������������.
	// TODO: ���� ����������� ����������� � ��������� 
	// ��������� ����� ������������ ���������� � ������
	int interpolationType;

	std::vector<std::vector<std::vector<double>>> EA_par;

	// �����
	std::vector<double> Energies;
	std::vector<double> Multiplicities;
};

enum particle_type { neutron = 0, proton, deutron, triton, alpha, recoils, gamma };

class mcEndfProduct
{
public:

	mcEndfProduct();

	~mcEndfProduct();

	//To do: enum n, gamma, proton
	particle_type product_type;

	int AWP, ZAP;
	
	std::vector<mcEndfEANuclearCrossSectionTable*> NuclearMultiplicity;

	// ������-������� ������� � ����������� �� ������� ���������� ��������
	std::vector<mcEndfEANuclearCrossSectionTable*> EANuclearCrossSections;
	
	void Load(std::istream& is);
};

// Class that keeps cross sections for one atomic isotope
class mcEndfP
{
public:
	mcEndfP();

	~mcEndfP();

	// ��������� �������, ���������� ������� ��� � ������� ���.
	// ������������ ��� ���������� �������������.
	std::string ElementName;

	// ������� ����� ������� ��������� � ������� ������� � ����������� �� ������� �������� �������
	// TODO: �������� ��������� �������. �����������, �� ����� �� ��� ������� ��������� 
	// � ���������� � ����, ��� � ������� ������ ������������� � dE/dX.
	mcEndfCrossSectionTable TotalCrossSections;

	// ������� ������� ������� � ����������� �� ������� �������� �������
	mcEndfCrossSectionTable NuclearCrossSections;

	std::vector<mcEndfProduct*> Products;

	// �������� ������ ����� �������
	void Load(const char* fname, const char* ename);

	void Clear();

	void dumpTotalCrossections(std::ostream& os) const;
};
