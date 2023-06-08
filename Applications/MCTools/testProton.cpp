#include <iostream>
#include "mcEndfP.h"
#include <fstream>
#include <filesystem>
#include "mcRng.h"
#include "mcScoreTest.h"


using namespace std;

void testproton() {
	const char* element = "Be009";
	std::string fname("../data/ENDFP/p-");
	fname += element;
	fname += ".dat";

	mcRng rng1, rng2;
	rng1.init(21, 48);
	mcEndfP elementData;
	elementData.Load(fname.c_str(), element);
	//elementData.dumpTotalCrossections(std::cout);

	//Testing incedent energy of proton, eV
	double kE = 2.2 * 1000000;
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
	mcRng rng;
	rng.init(33, 97);
	mcSpectrum A;
	A.init_inEn(kE);
	A.init_ptype(particle_type::proton);
	double** pars;
	double mu;
	double sigma = elementData.NuclearCrossSections.get_sigma(kE) * pow(10, -24);
	double ro = 7.874; // g/cm^3
	double const Na = 6.022 * pow(10, 23);
	double macro = ro * Na / 56 * sigma;
	double t = 7.5; // cm
	double probability = 1 - exp(-macro * t);
	for (int i = 0; i < 100000; i++)
	{
		if (probability > rng.rnd())
		{
			//n
			pars = elementData.EmittedNeutrons[0]->EANuclearCrossSections[0]->playpar(rng, kE, elementData.EmittedNeutrons[0]->LAW);
			mu = elementData.EmittedNeutrons[0]->EANuclearCrossSections[0]->playmu(kE, elementData.EmittedNeutrons[0]->LAW, pars, elementData.EmittedNeutrons[0]->product_type, rng);
			A.fill(pars[0][0], elementData.Products[0]->EANuclearCrossSections[0]->getMulti(kE), mu, elementData.Products[0]->product_type, rng);
			//p
			pars = elementData.Products[1]->EANuclearCrossSections[0]->playpar(rng, kE, elementData.Products[1]->LAW);
			mu = elementData.Products[1]->EANuclearCrossSections[0]->playmu(kE, elementData.Products[1]->LAW, pars, elementData.Products[1]->product_type, rng);
			A.fill(pars[0][0], elementData.Products[1]->EANuclearCrossSections[0]->getMulti(kE), mu, elementData.Products[1]->product_type, rng);
			//gamma
			pars = elementData.Products[elementData.Products.size() - 1]->EANuclearCrossSections[0]->playpar(rng, kE, elementData.Products[elementData.Products.size() - 1]->LAW);
			mu = elementData.Products[elementData.Products.size() - 1]->EANuclearCrossSections[0]->playmu(kE, elementData.Products[elementData.Products.size() - 1]->LAW, pars, elementData.Products[elementData.Products.size() - 1]->product_type, rng);
			A.fill(pars[0][0], elementData.Products[elementData.Products.size() - 1]->EANuclearCrossSections[0]->getMulti(kE), mu, elementData.Products[elementData.Products.size() - 1]->product_type, rng);
		}
	}
	A.dump(std::cout);
}