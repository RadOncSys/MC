#include "mcPhysicsNeutron.h"
#include "mcMediumNeutron.h"
#include "mcPhysicsCommon.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcDefs.h"
#include <float.h>

string CLEARFROMALPHA_(string x);

mcPhysicsNeutron::mcPhysicsNeutron(void)
{
}

mcPhysicsNeutron::~mcPhysicsNeutron(void)
{
}

bool mcPhysicsNeutron::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	if (p->ke <= ((const mcMediumNeutron&)med).transCutoff_neutron)
	{
		edep = p->ke;
		DiscardParticle(p);
		return true;
	}
	else
		return false;
}

double mcPhysicsNeutron::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	const mcMediumNeutron& m = (const mcMediumNeutron&)med;
	double logKE = ke;
	int iLogKE = int(ke);
	double sigma = (m.sigma0_neutro[iLogKE] + (logKE - iLogKE) * m.sigma1_neutro[iLogKE]) * dens;
	return (sigma > 0.0) ? 1 / sigma : DBL_MAX;
}

double mcPhysicsNeutron::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	p->p += p->u * step;
	p->dnear -= step;
	return 0;
}

double mcPhysicsNeutron::DoInterruction(mcParticle* p, const mcMedium* med) const
{

	mcRng& rng = p->thread_->rng();
	const mcMediumNeutron* m = (const mcMediumNeutron*)med;
	double logKE = p->ke;//log(ke);
	int iLogKE = int(p->ke);//int(m.iLogKE0_proto + logKE * m.iLogKE1_proto);
	double microsigma_total = (m->sigma0_neutro[iLogKE] + (logKE - iLogKE) * m->sigma1_neutro[iLogKE]) / m->density_ / NAVOGADRO * m->AtomicWeight();
	vector<double> sigmaratio;
	vector<double> probability;
	double psum = 0;
	for (int i = 0; i < m->elements_.size(); i++)
	{
		sigmaratio.push_back(m->elements_[i].partsByNumber * (m->Nmicrosigmaforelement(ROUND(m->elements_[i].atomicMass), ROUND(m->elements_[i].atomicNumber), p->ke, 0)) / microsigma_total);

		psum += sigmaratio[i];
	}
	for (int i = 0; i < m->elements_.size(); i++)
	{
		if (psum != 0)
			probability.push_back(sigmaratio[i] / psum);
		else break;
	}
	for (int i = 1; i < probability.size(); i++)
		probability[i] += probability[i - 1];
	double random = rng.rnd();
	int nucID = 0;
	for (nucID = 0; nucID < m->elements_.size(); nucID++)
	{
		if (random <= probability[nucID])
		{
			if (probability[nucID] == 0)
				continue;
			else break;
		}
	}
	//Теперь реакция осуществляется на nucID-ом ядре
	int endfID = 0;
	int A = ROUND(m->elements_[nucID].atomicMass);
	int Z = ROUND(m->elements_[nucID].atomicNumber);
	string elName = to_string(Z);
	if (A < 10)
		elName += "00" + to_string(ROUND(m->elements_[nucID].atomicMass));
	else if (A < 100)
		elName += "0" + to_string(A);
	else elName += to_string(A);
	for (endfID = 0; endfID < m->ENDFdata->size(); endfID++)
	{
		if (CLEARFROMALPHA_(m->ENDFdata->at(endfID)->ElementName) == elName)
		{
			break;
		}
	}
	//Теперь реакция осуществляется на endfID-ом ядре

	double El = 0, Inel = 0;
	std::vector<double> InelLVL;
	El = m->Nmicrosigmaforelement(A, Z, p->ke, 1);
	Inel = m->Nmicrosigmaforelement(A, Z, p->ke, 2);
	double Tot = El + Inel;
	El /= Tot;
	if (rng.rnd() < El)
		DoElastic(rng, endfID, p, m);
	else
	{
		for (int i = 0; i < m->ENDFdata->at(endfID)->nInelasticCS.size(); i++)
			InelLVL.push_back(m->Nmicrosigmaforelement(A, Z, p->ke, m->ENDFdata->at(endfID)->nInelasticCS[i]->MT));

		for (int i = 1; i < InelLVL.size(); i++)
		{
			InelLVL[i] += InelLVL[i - 1];
		}
		for (int i = 0; i < InelLVL.size(); i++)
		{
			InelLVL[i] /= InelLVL[InelLVL.size() - 1];
		}
		int LVLid = 0;
		double ksi1 = rng.rnd();
		for (LVLid = 0; LVLid < InelLVL.size(); LVLid++)
			if (InelLVL[LVLid] > ksi1)
				break;
	}
	double edep = p->ke / 2;
	p->ke = 0.0;
	return edep * p->weight;
}

void mcPhysicsNeutron::DoElastic(mcRng& rng, int endfID, mcParticle* p, const mcMediumNeutron* pmed)
{
	int i = 0;
	bool isLegendre = false;
	for (i = 0; i < pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies.size(); i++)
	{
		if (p->ke * 1000000 < pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies[i])
			break;
	}
	if (i < pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies.size())
		isLegendre = true;
	else for (i = 0; i < pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies.size(); i++)
	{
		if (p->ke * 1000000 < pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies[i])
			break;
	}
	if (isLegendre)
	{
		double cosCM = pmed->ENDFdata->at(endfID)->nElasticAngular.LegendreScat(i, rng);
	}
}

void mcPhysicsNeutron::DoInelastic()
{
}