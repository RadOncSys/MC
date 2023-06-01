#include <iostream>
#include "mcEndfP.h"
#include <fstream>
#include <filesystem>
#include "mcRng.h"
#include "mcScoreTest.h"


using namespace std;

void testproton() {
	const char* element = "Fe056";
	std::string fname("../data/ENDFP/p-");
	fname += element;
	fname += ".dat";

	mcRng rng1, rng2;
	mcEndfP elementData;
	elementData.Load(fname.c_str(), element);
	//elementData.dumpTotalCrossections(std::cout);

	//Testing incedent energy of proton, eV
	double kE = 100 * 1000000;
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
	for (int i = 0; i < 1000; i++)
	{
		//n
		pars = elementData.Products[0]->EANuclearCrossSections[0]->playpar(rng, kE);
		mu = elementData.Products[0]->EANuclearCrossSections[0]->playmu(kE, pars, elementData.Products[0]->product_type, rng);
		A.fill(pars[0][0], elementData.Products[0]->EANuclearCrossSections[0]->getMulti(kE), elementData.Products[0]->product_type);
		//p
		pars = elementData.Products[1]->EANuclearCrossSections[0]->playpar(rng, kE);
		mu = elementData.Products[1]->EANuclearCrossSections[0]->playmu(kE, pars, elementData.Products[1]->product_type, rng);
		A.fill(pars[0][0], elementData.Products[1]->EANuclearCrossSections[0]->getMulti(kE), elementData.Products[1]->product_type);
		//gamma
		pars = elementData.Products[72]->EANuclearCrossSections[0]->playpar(rng, kE);
		mu = elementData.Products[72]->EANuclearCrossSections[0]->playmu(kE, pars, elementData.Products[72]->product_type, rng);
		A.fill(pars[0][0], elementData.Products[72]->EANuclearCrossSections[0]->getMulti(kE), elementData.Products[72]->product_type);
	}
	A.dump(std::cout);
}