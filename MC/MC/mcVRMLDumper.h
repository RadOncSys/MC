// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <iostream>
using namespace std;

class mcVRMLDumper
{
public:
	mcVRMLDumper() { }
	~mcVRMLDumper(void) { }

	// Заголовок VRML файла. Выводится один раз в качестве первой записи файла.
	static void dumpHead(ostream& os);

	// оси координат мировой системы, a - длина осей (50 см по умолчанию)
	static void dumpWorldAxis(ostream& os, double a = 50.);
};
