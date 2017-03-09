#include "mcPhysicsCommon.h"
#include "mcDefs.h"
#include "mcRng.h"
#include <math.h>

//=================================
// Standard Physical functions
//---------------------------------

// ���������� ������� ������ ���� ��� ���������� ������� c targetMass � 
// ���������� � projectileMass � ������ �������� totalProjectileEnergy
// �� � MeV
// rpp p.321, eq.38.3
double ECenterOfMass(double targetMass, double projectileMass, double totalProjectileEnergy) {
	double a = SQUARE(targetMass) + SQUARE(projectileMass) + 2 * targetMass*totalProjectileEnergy;
	a = sqrt(a);
	return a;
}

//=================================
// �������� �������������� �������
//---------------------------------

// �������������� ������ gamma=sqrt(1/(1-beta^2))
// tEnergy - ������ ������� �������, Mass - � �����
inline double gamma(double Mass, double tEnergy)
{
	return tEnergy / Mass;
}

// �������� � �������� �������� ����� beta_squared=(velocity/light velocity)^2
// tEnergy - ������ ������� �������, Mass - � �����
double betasq(double Mass, double tEnergy) // 
{
	double	g_1 = Mass / tEnergy;
	return	1 - SQUARE(g_1);
}

// (beta*gamma) squared = P/M
// tEnergy - ������ ������� �������, Mass - � �����
double betagammasq(double Mass, double tEnergy)
{
	double	g = gamma(Mass, tEnergy);
	return	SQUARE(g) - 1;
}

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
int GaussStandardRnd_by_Marsaglia(mcRng& rng, double& r1, double& r2)
{
	double x, y, s2;
	do {
		x = 2 * rng.rnd() - 1;
		y = 2 * rng.rnd() - 1;
		s2 = SQUARE(x) + SQUARE(y);
	} while (s2 >= 1);
	r1 = sqrt(-2 * log(s2) / s2);
	r2 = y*r1;
	r1 = x*r1;
	return 0;
}