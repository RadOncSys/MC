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

	// Импорт данных протонов из множества файлов
	// fname - файл, содержащий заголовки сред, соглсованные с PEGS4
	// pstardir - файлы с тормозными способностями отдельных атомов по базе данных PSTAR
	// icrudir63 - файлы сечений ядерных реакций из протокола ICRU63
	void initProtonFromFiles(const string& fname, const string& nucleardir);
	void initProtonDeDxFromStream(istream&);
	void initProtonCSFromVector(std::shared_ptr<std::vector<std::shared_ptr<mcEndfP>>> dbData);

	void initNeutronFromStream(istream&);
	void initNeutronFromFiles(const string& fname, const string& nuclearDir);
	void initNeutronCSFromVector(std::shared_ptr<std::vector<std::shared_ptr<mcEndfN>>> dbData);

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

struct Mendeleev
{
	vector<bool> isNecessary;
	vector<bool> isLoad;

	void init();
};