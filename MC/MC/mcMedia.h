// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>
#include <memory>
using namespace std;

class mcMedium;
class mcMediumXE;
class mcMediumProton;
class mcMediumNeutron;
class mcPhysics;
class mcEndfP;
class mcEndfN;

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
	void initProtonFromFiles(const string& fname, const string& nucleardir);
	void initProtonDeDxFromStream(istream&);
	void initProtonCSFromVector(std::shared_ptr<std::vector<std::shared_ptr<mcEndfP>>> dbData);

	void initNeutronFromStream(istream&);
	void initNeutronFromFiles(const string& fname, const string& nuclearDir);
	void initNeutronCSFromVector(std::shared_ptr<std::vector<std::shared_ptr<mcEndfN>>> dbData);

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

struct Mendeleev
{
	vector<bool> isNecessary;
	vector<bool> isLoad;

	void init();
};