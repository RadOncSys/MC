#include <iostream>
#include "mcEndfP.h"
#include <fstream>
#include <filesystem>
#include "mcRng.h"
#include "mcScoreTest.h"

using namespace std;

void testEndf() 
{
	const char* element = "O016";
	std::string fname("../data/ENDFP/p-");
	fname += element;
	fname += ".tendl";

	mcEndfP elementData;
	elementData.Load(fname.c_str(), element);
	elementData.dumpTotalCrossections(std::cout);
}
