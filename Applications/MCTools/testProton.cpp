#include <iostream>
#include "mcEndfP.h"
#include <fstream>
#include <filesystem>
#include "mcRng.h"
#include "mcScoreTest.h"
#include "mcTransport.h"


using namespace std;

void testproton() {
	const char* element = "O016";
	std::string fname("../data/ENDFP/p-");
	fname += element;
	fname += ".dat";

	mcRng rng1, rng2;
	rng1.init(21, 48);
	mcEndfP elementData;
	elementData.Load(fname.c_str(), element);
	//elementData.dumpTotalCrossections(std::cout);
	//Testing incedent energy of proton, eV
	double kE = 72 * 1000000;
	//double I1 = elementData.Products[0]->EANuclearCrossSections[0]->integrate_f0(rng1, kE);
	double** pars;
	double mu;
	mcRng rng;
	rng.init(33, 97);
	for (int i = 0; i < 100; i++)
	{
		pars = elementData.Products[29]->EANuclearCrossSections[0]->playpar(rng, kE, elementData.Products[29]->LAW);
		cout << pars[0][0] << "  " << pars[1][0] << "  ";// << pars[2][0] << endl;
		mu = elementData.Products[29]->EANuclearCrossSections[0]->playmu(kE, elementData.Products[29]->LAW, pars, elementData.Products[29]->product_type, rng);
		cout << mu;
		cout << endl;
	}
	cout << endl;
	/*rng1.init(55, 97);
	rng2.init(40, 97);*/

	//cout << "Getmulti = " << elementData.Products[0]->EANuclearCrossSections[0]->getMulti(kE) << endl;
	/*for (int p_i = 0; p_i < 10; p_i++)
	{
		vector<vector<double>> Multi_;
		Multi_.resize(6);
		for (int i = 0; i < Multi_.size(); i++)
		{
			Multi_[i].resize(3);
			Multi_[i][0] = 0;
		}
		double Norm = 0;
		int c = 0;
		for (int i = 0; i < elementData.Products.size(); i++)
		{
			int a = elementData.Products[i]->product_type;
			if (a <= 4 || a == 6)
			{
				Multi_[c][1] = elementData.Products[i]->EANuclearCrossSections[0]->getMulti(kE);
				for (int j = c; j < Multi_.size(); j++)
				{
					Multi_[j][0] += elementData.Products[i]->EANuclearCrossSections[0]->getMulti(kE);
					Multi_[j][2] = a;
				}
				Norm += elementData.Products[i]->EANuclearCrossSections[0]->getMulti(kE);;
				c++;
			}
		}
		double r1 = rng1.rnd();
		double r2 = rng2.rnd();
		for (int i = 0; i < Multi_.size(); i++)
		{
			Multi_[i][0] /= Norm;
		}*/
		//for (int ii = 0; ii < Multi_.size(); ii++)
		//{
		//	if (r1 < Multi_[ii][0])
		//	{
		//		cout << "#" << p_i + 1 << endl;
		//		if (Multi_[ii][1] > 1 && Multi_[ii][1] < 2) 
		//			if (r2 > Multi_[ii][1] - 1)
		//				cout << "One ";		
		//			else cout << "Two ";
		//		if (Multi_[ii][1] > 2)
		//			if (r2 > Multi_[ii][1] - 2)
		//				cout << "Two ";
		//			else cout << "Three ";
		//		cout << typeof(Multi_[ii][2]) << " is out." << endl;
		//		break;
		//	}
		//}
	//}
				//mcRng rng;
				//rng.init(33, 97);
				//mcSpectrum A;
				//A.init_inEn(kE);
				//A.init_ptype(particle_type::proton);
				//double** pars;
				//double mu;
				//double ro = 1.00; // g/cm^3
				//double lambda = elementData.NuclearCrossSections.get_lambda(kE, ro, 18);
				//double Sigma = 1 / lambda;
				//vector<double> inteructDist;
				//double meantointeruct = 0;
				//for (int i = 0; i < 1000000; i++)
				//{
				//	inteructDist.push_back(mcTransport::HowManyMFPs(rng) * lambda);				//ÇÀÄÀ×À ÂÀËÈÄÀÖÈÈ #1 äëÿ Âîäû - ðåøåíà
				//	meantointeruct += inteructDist[i];
				//}
				//meantointeruct /= inteructDist.size();
				//double t = 7.5; // cm
				//double probability = 1 - exp(-t / lambda);
	//for (int i = 0; i < 100000; i++)
	//{
	//	if (probability > rng.rnd())
	//	{
	//		//n
	//		pars = elementData.EmittedNeutrons[0]->EANuclearCrossSections[0]->playpar(rng, kE, elementData.EmittedNeutrons[0]->LAW);
	//		mu = elementData.EmittedNeutrons[0]->EANuclearCrossSections[0]->playmu(kE, elementData.EmittedNeutrons[0]->LAW, pars, elementData.EmittedNeutrons[0]->product_type, rng);
	//		A.fill(pars[0][0], elementData.Products[0]->EANuclearCrossSections[0]->getMulti(kE), mu, elementData.Products[0]->product_type, rng);
	//		//p
	//		pars = elementData.Products[1]->EANuclearCrossSections[0]->playpar(rng, kE, elementData.Products[1]->LAW);
	//		mu = elementData.Products[1]->EANuclearCrossSections[0]->playmu(kE, elementData.Products[1]->LAW, pars, elementData.Products[1]->product_type, rng);
	//		A.fill(pars[0][0], elementData.Products[1]->EANuclearCrossSections[0]->getMulti(kE), mu, elementData.Products[1]->product_type, rng);
	//		//gamma
	//		pars = elementData.Products[elementData.Products.size() - 1]->EANuclearCrossSections[0]->playpar(rng, kE, elementData.Products[elementData.Products.size() - 1]->LAW);
	//		mu = elementData.Products[elementData.Products.size() - 1]->EANuclearCrossSections[0]->playmu(kE, elementData.Products[elementData.Products.size() - 1]->LAW, pars, elementData.Products[elementData.Products.size() - 1]->product_type, rng);
	//		A.fill(pars[0][0], elementData.Products[elementData.Products.size() - 1]->EANuclearCrossSections[0]->getMulti(kE), mu, elementData.Products[elementData.Products.size() - 1]->product_type, rng);
	//	}
	//}
	//A.dump(std::cout);
}