#pragma once
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include "mcRng.h"
#include "mcEndfP.h"

using namespace std;

class mcSpectrum
{
public:
	//����� ����������� ���� ������ �������
	void fill(double Eout, double weight, int p_type);

	//����� �������
	void dump(std::ostream& os) const;

	//���������������� ������ (E*n �� E):
	//ESpectrum[0] - n
	//Espectrum[1] - p
	//Espectrum[2] - gamma
	//ESpectrum[i][j] - E * n for corresponding particle, j is from 0 to 150 (MeV)
	vector<vector<double>> ESpectrum;

	vector<double> Energies;
	
	bool energy_resized = false;

	//������� ���������� ������
	double inEn;

	void init_inEn(double kE) {
		inEn = kE;
	}

	//��� ���������� ������
	particle_type inP;

	void init_ptype(particle_type ptype) {
		inP = ptype;
	}

};