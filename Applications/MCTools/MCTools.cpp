// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
// Auxilary project to make complex tests of MC code
//---------------------------------------------------------------------------
#include <iostream>
#include "mcEndfP.h"

void testproton();

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::wcout << L"Usage:" << std::endl;
		std::wcout << argv[0] << L" mode" << std::endl << std::endl;
		std::wcout << L"modes:" << std::endl;
		std::wcout << "   1 - " << L"dump example nuclear reaction crossections for protons" << std::endl;
		std::wcout << std::endl;
		return -1;
	}

	try
	{
		if (strcmp(argv[1], "1") == 0)
		{
			const char* element = "O016";
			std::string fname("../data/ENDFP/p-");
			fname += element;
			fname += ".tendl";

			mcEndfP elementData;
			elementData.Load(fname.c_str(), element);
			elementData.dumpTotalCrossections(std::cout);
		}
		else if (strcmp(argv[1], "3") == 0)
		{
			testproton();
		}
	}
	catch (std::exception& e)
	{
		std::cout << "Tool error: " << e.what() << std::endl;
	}
	catch (...)
	{
		std::cout << "Unknown tool error ..." << std::endl;
	}

}
