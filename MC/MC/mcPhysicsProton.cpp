#include "mcPhysicsProton.h"
#include "mcMediumProton.h"
#include "mcPhysicsCommon.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcDefs.h"
#include <float.h>

mcPhysicsProton::mcPhysicsProton(void)
{
}

mcPhysicsProton::~mcPhysicsProton(void)
{
}

bool mcPhysicsProton::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	if (p->ke <= ((const mcMediumProton&)med).transCutoff_proto)
	{
		edep = p->ke;
		DiscardParticle(p);
		return true;
	}
	else
		return false;
}

double mcPhysicsProton::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	// ������ ��������� ���������� ������� �������
	const mcMediumProton& m = (const mcMediumProton&)med;
	double logKE = ke;//log(ke);
	int iLogKE = int(ke);//int(m.iLogKE0_proto + logKE * m.iLogKE1_proto);
	double sigma = (m.sigma0_proto[iLogKE] + logKE * m.sigma1_proto[iLogKE]) * dens;
	return (sigma > 0.0) ? 1 / sigma : DBL_MAX;
}

double mcPhysicsProton::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	const mcMediumProton& m = (const mcMediumProton&)med;
	double e_dep = 0;
	// VK
	// ���������� ������ ����, �� ������� ���������� ��������� � ������� ��� � step
	// ���������� ������ ������� �� ���� ���� 
	// � ������ ����������
	// � ������� �� return e_dep 
	// �������� ��������� �������
	// ��������� ��������� - �������� �����������

	// 1. ���������� ������� ������ �������
	// ���� ��������, ��� ������ �� ������ ���� ������ �� ������ 1% �������
	e_dep = 0.01 * p->ke;	// ����� ������������ �������� ������ ������� �� ����
	// �������� �� ����� ������� ����� ����
	// ���� dE/dx ��� ������� ������� ����� ����: kE-dE/2.0
	//double mE = kE;		// �������� ������ (�����)
	double mE = p->ke - e_dep / 2.0;  // ��������� ���� 2008
	// ������ ��������� ���������� ������� �������
	double logKE = mE; //p->ke;//log(p->ke);
	int iLogKE = int(logKE);//(int) (m.iLogKE0_proto + logKE * m.iLogKE1_proto);
	double dedx = p->regDensityRatio * (m.dedx0_proto[iLogKE] + logKE * m.dedx1_proto[iLogKE]);
	// ����� �������� ���� (�� ����������� � �������� ��� ��������� ��. ������� ��� ��������
	step = MIN(step, e_dep / dedx);
	e_dep = step * dedx;	// ���� ��������� ���

	// �������� �� ��, ��� ������ ������ ����� ���� �� ������
	// 2. ��������� ����������
	double rel = (p->ke + PMASS) * PMASS_1;
	double sig2 = p->regDensityRatio * m.dEdxStragglingGaussVarianceConstPart_ * step * (1 + SQUARE(rel));

	double r1, r2;
	GaussStandardRnd_by_Marsaglia(p->thread_->rng(), r1, r2); // �������� ������������� �� ������ ��������
	e_dep = e_dep - sqrt(sig2) * r1;
	e_dep = (e_dep < 0) ? 0 : (e_dep > p->ke) ? p->ke : e_dep;

	// 3. �������������� �������
	p->p += p->u * step;
	//p->ke -= e_dep;

	// 4. ����������� ����� ����������� � ������������� ���
	// ���������� ������������ ���������
	double tE = p->ke + PMASS; // ������ �������
	// �� ������ �� ���� p->regDensityRatio, �� � Nova ����� ������ ��� ������
	double l = p->regDensityRatio * step / m.radLength;
	//double a1	= 1+0.038*log(l); // ��������� � �������
	//double a2	= betasq(PMASS,tE)*tE; // ����������� b*c*p=E*b^2
	//double a1_2 = a1/a2;
	double a1_2 = (1 + 0.038 * log(l)) / (betasq(PMASS, tE) * tE);
	double th0 = sqrt(184.96 * l * SQUARE(a1_2)); // 184.96=13.6^2	
	// ��-���� �� �� ����� ����, � � EGS ������������ �� ������. ����� ������ �������
	// ���������� theta � ���������� ��
	// ��. rpp-2006-book.pdf 27.3 p.262 eq.27.18,27.19
	// � ������ ���� ���� ����������� �� 2*pi
	// ���� �� ����������� � ���, ��� � Nova � ��� ��� �������� ��� p+ (�����!)
	// ��� ��������� �����, ��� �����
	double theta = fabs(r2 * th0); //eq.27.19
	//�������� ����������
	//double sqrtOf12inverse = 0.28867513459481288225457439025098;	
	//double deflection	= r2*th0*t*sqrtOf12inverse+theta*t*0.5;	//eq.27.18
	double phi = p->thread_->rng().rnd() * TWOPI;
	ChangeDirection(cos(theta), sin(theta), cos(phi), sin(phi), p->u);

	// p->mfps = 0; // VK ��� �����? ���� �� �� ����������� ������ ��� mfp, � ���������� ���, �� ��� ������ ������� �� ����� 

	p->ke -= e_dep;
	return e_dep;
}

double mcPhysicsProton::DoInterruction(mcParticle* p, const mcMedium* med) const
{
	// TODO: ������� ��������� ����������� � ����������� ����������������.
	// 1. ��������� ����� ���������� �������������� �� ������ ��������.
	//    �������� �� ����� ������������ ���������, ��� ��� ���������� 2 �������
	//    ���������� ������� ������ � ������� ��������.
	// 2. � ���������� ������� �������������� ������ ������ �����������.
	//    � ���� �� ������������ ������� �������, �� ����-�� ������������ � 
	//    ��������� ������������� ���������� ��������, ������ � ����������� �� ��
	//    ������� �� ��������������.

	// ��� ���������� �������������� ���������� ��������� �������.
	// � ������� �������� 1/3 ������������ ������� (� ���� �� �������� ???).
	// 1/3 ����������� � �����.
	// 1/3 ��������� �� ������� ������� ��������.

	p->ke /= 3.0;
	//cout << ".";

	// ���������� �������, ������������ � �����.
	return 2 * p->ke * p->weight;
}
