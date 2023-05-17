#include <iostream>
#include "mcEndfP.h"

void testproton() {
	const char* element = "Fe056";
	std::string fname("../data/ENDFP/p-");
	fname += element;
	fname += ".dat";

	mcEndfP elementData;
	elementData.Load(fname.c_str(), element);
	elementData.dumpTotalCrossections(std::cout);
}