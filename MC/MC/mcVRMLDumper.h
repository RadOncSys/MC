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

	// ��������� VRML �����. ��������� ���� ��� � �������� ������ ������ �����.
	static void dumpHead(ostream& os);

	// ��� ��������� ������� �������, a - ����� ���� (50 �� �� ���������)
	static void dumpWorldAxis(ostream& os, double a = 50.);
};
