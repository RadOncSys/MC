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
	double edep = 0;
	/*
	mcRng& rng = p->thread_->rng();
	const mcMediumNeutron* m = (const mcMediumNeutron*)med;
	double logKE = p->ke;//log(ke);
	int iLogKE = int(p->ke);//int(m.iLogKE0_proto + logKE * m.iLogKE1_proto);
	double microsigma_total = (m->sigma0_neutro[iLogKE] + (logKE - iLogKE) * m->sigma1_neutro[iLogKE]) / m->density_ / NAVOGADRO * m->AtomicWeight();
	if (microsigma_total <= 0)
		return 0;
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
	{
		double ke_before = p->ke;
		DoElastic(rng, endfID, p, m, A);
		edep = ke_before - p->ke;
	}
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
		if (LVLid < m->ENDFdata->at(endfID)->nInelasticCS.size() - 1 || m->ENDFdata->at(endfID)->nInelasticCS.back()->MT != 91)
		{
			double ke_before = p->ke;
			DoInelastic(rng, endfID, LVLid, p, m, A);
		}
		else
		{
			double ke_before = p->ke;
			DoInelasticCont(rng, endfID, LVLid, p, m);
			edep = ke_before - p->ke;
		}
	}
	*/
	return edep * p->weight;
}

void mcPhysicsNeutron::DoElastic(mcRng& rng, int endfID, mcParticle* p, const mcMediumNeutron* pmed, int A)
{
	/*
	int i = 0;
	bool isLegendre = false;
	double cosCM = 0;
	if (pmed->ENDFdata->at(endfID)->nElasticAngular.LI == 1)
	{
		cosCM = rng.rnd();
		if (rng.rnd() > 0.5)
			cosCM *= -1;
	}
	else
	{
		for (i = 0; i < pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies.size(); i++)
		{
			if (p->ke * 1000000 < pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies[i])
				break;
		}
		if (i < pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies.size() || i < pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies.size() == 0)
			isLegendre = true;
		else for (i = 0; i < pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies.size(); i++)
		{
			if (p->ke * 1000000 < pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies[i])
				break;
		}
		if (isLegendre)
		{
			if (i == pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies.size())
				i--;
			if (pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies[i] - p->ke > p->ke - pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies[i - 1] && i != 0)
				i--;
			cosCM = pmed->ENDFdata->at(endfID)->nElasticAngular.LegendreScat(i, rng);
		}
		else
		{
			if (i == pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies.size())
				i--;
			if (pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies[i] - p->ke > p->ke - pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies[i - 1] && i != 0)
				i--;
			cosCM = pmed->ENDFdata->at(endfID)->nElasticAngular.TableScat(i, rng);
		}
	}
	double ke_ = p->ke * (A * A + 2 * A * cosCM + 1) / (1 + 2 * A + A * A);
	double cosLS = (A + 1) / 2 * sqrt(ke_ / p->ke) - (A - 1) / 2 * sqrt(p->ke / ke_);
	double sinphi = sin(2 * PI * rng.rnd());
	double cosphi = cos(2 * PI * rng.rnd());
	ChangeDirection(cosLS, sin(acos(cosLS)), cosphi, sinphi, p->u);
	p->ke = ke_;
	*/
}

void mcPhysicsNeutron::DoInelastic(mcRng& rng, int endfID, int LVLid, mcParticle* p, const mcMediumNeutron* pmed, int A)
{
	/*
	int i = 0;
	bool isLegendre = false;
	double cosCM = 0;
	if (pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->LI == 1)
	{
		cosCM = rng.rnd();
		if (rng.rnd() > 0.5)
			cosCM *= -1;
	}
	else
	{
		for (i = 0; i < pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->LEnergies.size(); i++)
		{
			if (p->ke * 1000000 < pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->LEnergies[i])
				break;
		}
		if (i < pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->LEnergies.size() || i < pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies.size() == 0)
			isLegendre = true;
		else for (i = 0; i < pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->TEnergies.size(); i++)
		{
			if (p->ke * 1000000 < pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->TEnergies[i])
				break;
		}
		if (isLegendre)
		{
			if (i == pmed->ENDFdata->at(endfID)->nElasticAngular.LEnergies.size())
				i--;
			if (pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->LEnergies[i] - p->ke > p->ke - pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->LEnergies[i - 1] && i != 0)
				i--;
			cosCM = pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->LegendreScat(i, rng);
		}
		else
		{
			if (i == pmed->ENDFdata->at(endfID)->nElasticAngular.TEnergies.size())
				i--;
			if (pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->TEnergies[i] - p->ke > p->ke - pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->TEnergies[i - 1] && i != 0)
				i--;
			cosCM = pmed->ENDFdata->at(endfID)->inelasticLevelsAng[LVLid]->TableScat(i, rng);
		}
	}
	double Q = pmed->ENDFdata->at(endfID)->nInelasticCS[LVLid]->Q;
	double alpha = (A * A - 2 * A + 1) / (A * A + 2 * A + 1);
	double ke_ = p->ke / 2 * (1 + alpha - 2 * A / (1 + A) * Q / p->ke + (1 - alpha) * cosCM * sqrt(1 - (A + 1) / A * Q / p->ke));
	double cosLS = (1 + A * cosCM) / sqrt(1 + 2 * A * cosCM + A * A);
	double sinphi = sin(2 * PI * rng.rnd());
	double cosphi = cos(2 * PI * rng.rnd());
	ChangeDirection(cosLS, sin(acos(cosLS)), cosphi, sinphi, p->u);
	p->ke = ke_;
	*/
}

void mcPhysicsNeutron::DoInelasticCont(mcRng& rng, int endfID, int LVLid, mcParticle* p, const mcMediumNeutron* pmed)
{
	/*
	double primary_ke = p->ke;
	int Nquantity = pmed->ENDFdata->at(endfID)->nInelasticContin[0]->EANuclearCrossSections[0]->playMulti(p->ke * 1000000, rng);
	if (Nquantity > 1)
		throw exception("Multi neutron during inelastic scattering!?");
	int Gquantity = pmed->ENDFdata->at(endfID)->nInelasticContin[2]->EANuclearCrossSections[0]->playMulti(p->ke * 1000000, rng);
	if (pmed->ENDFdata->at(endfID)->nInelasticContin[0]->LAW != 1)
		throw exception("This LAW doesn't exist in ENDF play-block.");
	for (int i = 0; i < Gquantity; i++)
	{
		mcParticle* pNewPhoton = DuplicateParticle(p);
		pNewPhoton->t = MCP_PHOTON;
		pNewPhoton->q = 0;
		int eoutID = 0, keIN = 0;
		pNewPhoton->ke = pmed->ENDFdata->at(endfID)->nInelasticContin[2]->EANuclearCrossSections[0]->playE(primary_ke, keIN, eoutID, rng);
		GoInRandomDirection(rng.rnd(), rng.rnd(), pNewPhoton->u);
		p->ke -= pNewPhoton->ke;
	}
	mcParticle* pNewNeutron = DuplicateParticle(p);
	pNewNeutron->t = MCP_NEUTRON;
	pNewNeutron->q = 0;
	int eoutID = 0, keIN = 0;
	double neutron_ke = pmed->ENDFdata->at(endfID)->nInelasticContin[0]->EANuclearCrossSections[0]->playE(primary_ke, keIN, eoutID, rng);
	getKallbachMannAngle(rng, endfID, pNewNeutron, pmed, keIN, eoutID);
	p->ke -= neutron_ke;
	pNewNeutron->ke = neutron_ke;
	*/
}

void mcPhysicsNeutron::getKallbachMannAngle(mcRng& rng, int endfID, mcParticle* p, const mcMediumNeutron* pmed, int keIN, int eoutID)
{
	/*
	double ke_ = pmed->ENDFdata->at(endfID)->nInelasticContin[0]->EANuclearCrossSections[0]->Energies[keIN]; //В первом приближении энергия без интерполяции
	double costheta = pmed->ENDFdata->at(endfID)->nInelasticContin[0]->EANuclearCrossSections[0]->playmu(ke_, pmed->ENDFdata->at(endfID)->nInelasticContin[0]->LAW, keIN, eoutID, 0, rng);
	double phi = 2 * PI * rng.rnd();
	double cosphi = cos(phi);
	double sinphi = sin(phi);
	double AWRa = 1, AWRA = pmed->ENDFdata->at(endfID)->nInelasticContin[0]->EANuclearCrossSections[0]->AWR_nucl;
	double AWRb = 0;
	double Eb = p->ke;
	ke_ /= 1000000;
	p->ke = Eb + AWRa * AWRb * ke_ / (AWRA + AWRa) / (AWRA + AWRa) + 2 * sqrt(AWRa * AWRb * ke_ * Eb) * costheta / (AWRA + AWRa);
	costheta = sqrt(Eb / p->ke) * costheta + sqrt(AWRa * AWRb * ke_ / p->ke) / (AWRA + AWRa);
	double sintheta = sin(acos(costheta));
	ChangeDirection(costheta, sintheta, cosphi, sinphi, p->u);
	return;
	*/
}