// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
// Classes to manage cross sections load for proton nuclear interractions.
// Current version supports ICRU-63 data.
//---------------------------------------------------------------------------
#pragma once
#include <string>
#include <vector>
#include <memory>

// ����� �������� ���������� ������ ��� �������� ������ ����
class mcCSNuclearForAngleSpectrum
{
public:
	// ���� � ��������, ��� �������� ����� ������ ������
	double Angle;

	// ������� ������� ��������� ������� ������� ���� �������
	std::vector<double> SourceBinUps;

	// ������������� ����� �������� �������
	std::vector<double> SourceSpectrum;

	void Clear();
};

// Crossections per istope per incident particle energy
class mcCSNuclearParticleEnergy
{
public:
	// ������������ ������� ��������� �������
	double Energy;

	// ��������� ������� ������� �������������� 
	// ������ ������� ��� ������������� ������� (mb)
	double TotalCrossSection;
	double ProtonCrossSection;
	double NeutronCrossSection;

	// ������ �������� ����������� ��������
	std::vector<mcCSNuclearForAngleSpectrum> ProtonAngles;

	// ������ �������� ����������� ���������
	std::vector<mcCSNuclearForAngleSpectrum> NeutronAngles;

	// Other particles
	// � ������ ������ �� ������������ � ������������� ���������� ��������� �������.
	// �.�., ������� ����������� ���������� ������� �������� � ��������� 
	// �������� ��������� ������ � ��������� ��� ������������ ������� � �����.

	void Clear();
};

// Class that keeps coss sections for one atomic isotope
class mcCSNuclear
{
public:
	// ��������� �������, ���������� ������� ��� � ������� ���.
	// ������������ ��� ���������� �������������.
	std::string ElementName;

	// ������� ������ ������� �������� ��������� �������
	std::vector<mcCSNuclearParticleEnergy> Energies;

	void Clear();
};
