#include <iostream>
#include "mcEndfNP.h"
#include <fstream>
#include <filesystem>
#include "mcRng.h"
#include "mcScoreTest.h"

using namespace std;

void testNeutrons()
{
	const char* element = "1H001";
	std::string fname("../data/ENDFN/n-");
	fname += element;
	fname += ".dat";

	mcEndfNP elementData;
	elementData.Load(fname.c_str(), element);
	elementData.dumpTotalCrossections(std::cout);
}
