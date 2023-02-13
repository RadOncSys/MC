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
class mcMediumProton;
class mcMediumNeutron;
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
	const mcMediumProton* getProtonMedium(short idx) const;
	const mcMediumProton* getNeutronMedium(short idx) const;

	void initXEFromStream(istream&);
	void initXEFromFile(const string& fname);

	// ������ ������ �������� �� ��������� ������
	// fname - ����, ���������� ��������� ����, ������������ � PEGS4
	// pstardir - ����� � ���������� ������������� ��������� ������ �� ���� ������ PSTAR
	// icrudir63 - ����� ������� ������� ������� �� ��������� ICRU63
	void initProtonFromFiles(const string& fname, const string& pstardir, const string& icru63dir);
	void initProtonDeDxFromStream(istream&);

	void initNeutronFromStream(istream&);
	void initNeutronFromFile(const string& fname);

	// ���������� ��������� ������� ���������� �������� ��� ������� ���������� ����
	const mcPhysics* getPhysics(int ptype) const;

	// ���������� ��������� ������� ���������� �������� ��� ������� ���������� ����
	const mcMedium* getMedium(int ptype, int idx) const;

	vector<mcMedium*> Media() { return xes_; }

protected:
	// ������ ���� ����
	vector<string> mnames_;

	// ��������� ���� ���������� �������, ����������, �������� � ���������
	vector<mcMedium*> xes_;
	vector<mcMedium*> protons_;
	vector<mcMedium*> neutrons_;

	vector<mcPhysics*> physics_;
};
