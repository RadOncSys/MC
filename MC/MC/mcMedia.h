// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>
using namespace std;

class mcMedium;
class mcMediumXE;
class mcPhysics;

class mcMedia
{
public:
	mcMedia(void);
	~mcMedia(void);

	// ���������� ���� �� ���������.
	// ���������� ������� ���������� ����� ������������ ���� �������
	// ������� ������������ �� ������ ��� ������, �������� �� ������� ����� ������.
	void addName(const char* mname);

	// ����� ����� �� �����, �������� ��� ������������� �����������
	short getMediumIdx(const char* mname) const;

	const mcMediumXE* getMediumXE(short idx) const;

	void initXEFromStream(istream&);
	void initXEFromFile(const string& fname);

	// ���������� ��������� ������� ���������� �������� ��� ������� ���������� ����
	const mcPhysics* getPhysics(int ptype) const;

	// ���������� ��������� ������� ���������� �������� ��� ������� ���������� ����
	const mcMedium* getMedium(int ptype, int idx) const;

	vector<mcMedium*> Media() { return xes_; }

protected:
	// ������ ���� ����
	vector<string> mnames_;

	// ��������� ���� ���������� ������� � ����������
	vector<mcMedium*> xes_;

	vector<mcPhysics*> physics_;
};
