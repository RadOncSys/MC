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

// Class that keeps coss sections for one atomic isotope
class mcEndfP
{
public:
	// ��������� �������, ���������� ������� ��� � ������� ���.
	// ������������ ��� ���������� �������������.
	std::string ElementName;

	// ������� ������� ������� � ����������� �� ������� �������� �������
	mcEndfCrossSectionTable CrossSections;

	// �������� ������ ����� �������
	void Load(const char* fname, const char* ename);

	void Clear();
};
