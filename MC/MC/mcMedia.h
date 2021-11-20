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

	// Добавление сред по названиям.
	// Подготовка сечений происходит после перечисления сред вызовом
	// функций иницилизации из потока или файлов, отдельно по группам типов частиц.
	void addName(const char* mname);

	// Поиск среды по имени, полезный для инициализации транспортов
	short getMediumIdx(const char* mname) const;

	const mcMediumXE* getMediumXE(short idx) const;
	const mcMediumProton* getProtonMedium(short idx) const;
	const mcMediumProton* getNeutronMedium(short idx) const;

	void initXEFromStream(istream&);
	void initXEFromFile(const string& fname);

	void initProtonFromStream(istream&);
	void initProtonFromFile(const string& fname);

	void iniNeutronFromStream(istream&);
	void iniNeutronFromFile(const string& fname);

	// Возвращает указатель объекта физических расчетов для частицы указанного типа
	const mcPhysics* getPhysics(int ptype) const;

	// Возвращает указатель объекта физических расчетов для частицы указанного типа
	const mcMedium* getMedium(int ptype, int idx) const;

	vector<mcMedium*> Media() { return xes_; }

protected:
	// Список имен сред
	vector<string> mnames_;

	// Параметры сред транспорта фотонов, электронов, протонов и нейтронов
	vector<mcMedium*> xes_;
	vector<mcMedium*> protons_;
	vector<mcMedium*> neutrons_;

	vector<mcPhysics*> physics_;
};
