// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcRng;

//=================================
// Standard Physical functions
//---------------------------------

// ���������� ������� ������ ���� ��� ���������� ������� c targetMass � 
// ���������� � projectileMass � ������ �������� totalProjectileEnergy
// �� � MeV
// rpp p.321, eq.38.3
double ECenterOfMass(double targetMass, double projectileMass, double totalProjectileEnergy);

//=================================
// �������� �������������� �������
//---------------------------------

// �������������� ������ gamma=sqrt(1/(1-beta^2))
// tEnergy - ������ ������� �������, Mass - � �����
inline double gamma(double Mass, double tEnergy);

// �������� � �������� �������� ����� beta_squared=(velocity/light velocity)^2
// tEnergy - ������ ������� �������, Mass - � �����
double betasq(double Mass, double tEnergy);

// (beta*gamma) squared = P/M
// tEnergy - ������ ������� �������, Mass - � �����
double betagammasq(double Mass, double tEnergy);

//=================================
// Standard Random Numbers functions
//---------------------------------
// �������� ������������� ��������� ��������
// �������� 2-� ��������� �������� r1 r2, ������������ ����������� ������
// (��������� �� ������� = 0 � ����������� ����������� 1)
// ������� Marsaglia
// ��� ������ ��� ������ ���������� ������
// ��. http://en.wikipedia.org/wiki/Marsaglia_polar_method
// ��������� ���������� ������������� ����� rng.rnd() ������ ������ [0,1] -?
int GaussStandardRnd_by_Marsaglia(mcRng& rng, double& r1, double& r2);
