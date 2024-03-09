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
		sigmaratio.push_back(m->elements_[i].partsByNumber * (m->Nmicrosigmaforelement(ROUND(m->elements_[i].atomicMass), ROUND(m->elements_[i].atomicNumber), p->ke)) / microsigma_total);

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

	double El = 0, Inel = 0;
	El = m->ENDFdata->at(endfID)->ElasticCrossSections.get_sigma(p->ke * 1000000);
	if (!m->ENDFdata->at(endfID)->InelasticCrossSections.isEmpty)
		Inel = m->ENDFdata->at(endfID)->InelasticCrossSections.get_sigma(p->ke * 1000000);
	double Tot = El + Inel;
	El /= Tot;
	if (rng.rnd() < El)
		DoElastic(rng, endfID, p, m);
	else DoInelastic();
	double edep = p->ke / 2;
	p->ke = 0.0;
	return edep * p->weight;
}

void mcPhysicsNeutron::DoElastic(mcRng& rng, int endfID, mcParticle* p, const mcMediumNeutron* pmed)
{
}

void mcPhysicsNeutron::DoInelastic()
{
}
