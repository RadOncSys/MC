// Radiation Oncology Monte Carlo open source project
//
// Author: [2023] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
// Classes to manage cross sections load for proton nuclear interractions
// from ENDF format
//---------------------------------------------------------------------------
#pragma once
#include <string>
#include <vector>
#include <memory>

// Структура стрроки ENDF файла
struct mcEndfRecord
{
	char c[6][11];
	//char data[66];
	char Z[2];
	char Stblt[2];
	char MF[2];
	char MT[3];
	char LineNumber[5];

	// Парсинг значени с плавающей точкой в формате ENDF
	// (где степнь указана нестандартно после знака +/-).
	static double ParseValue(const char* s, int n);
};

// Crossections per istope per incident particle energy
class mcEndfCrossSectionTable
{
public:
	void Load(std::istream& is);
	void dump(std::ostream& os) const;

	// Количество пар энергия падающей частицы / сечение
	int npoints;

	// Количество видов интерполяции
	// Временно предполагаем, что мы не столкнемся со 
	// множесвенностью интерполяций в пределах одной таблицы.
	// TODO: Проверить эту гипотезу (в коде стоит exception на этот случай)
	int ninterpolations;
	 
	// Тип интерполяции.
	// TODO: если обнаружится потребность в поддержке 
	// множества типов интерполяций переделать в массив
	int interpolationType;
	
	// Точки
	std::vector<double> Energies;
	std::vector<double> Values;
};

// Класс для чтения MF=6 MT=5 

class mcEndfEANuclearCrossSectionTable
{
public:
	void mLoad(std::istream& is);
	//void dump(std::ostream& os) const;

	void Load(std::istream& is);

	// Количество пар энергия падающей частицы / мультиплетность
	int n_energypoints;

	// Количество чего-то там
	int npoints_ang;

	// Количество видов интерполяции
	// Временно предполагаем, что мы не столкнемся со 
	// множесвенностью интерполяций в пределах одной таблицы.
	// TODO: Проверить эту гипотезу (в коде стоит exception на этот случай)
	int ninterpolations;

	// Тип интерполяции.
	// TODO: если обнаружится потребность в поддержке 
	// множества типов интерполяций переделать в массив
	int interpolationType;

	std::vector<std::vector<std::vector<double>>> EA_par;

	// Точки
	std::vector<double> Energies;
	std::vector<double> Multiplicities;
};

enum particle_type { neutron = 0, proton, deutron, triton, alpha, recoils, gamma };

class mcEndfProduct
{
public:

	mcEndfProduct();

	~mcEndfProduct();

	//To do: enum n, gamma, proton
	particle_type product_type;

	int AWP, ZAP;
	
	std::vector<mcEndfEANuclearCrossSectionTable*> NuclearMultiplicity;

	// Энерго-угловые сечения в зависимости от энергии налетающих протонов
	std::vector<mcEndfEANuclearCrossSectionTable*> EANuclearCrossSections;
	
	void Load(std::istream& is);
};

// Class that keeps cross sections for one atomic isotope
class mcEndfP
{
public:
	mcEndfP();

	~mcEndfP();

	// Назавание изотопа, включающее атомное имя и атомный вес.
	// Используется как уникальный идентификатор.
	std::string ElementName;

	// Сечения суммы упругих рассеяний и ядерных реакций в зависимости от энергии падающей частицы
	// TODO: Возможно временная таблица. Разобраться, не нужно ли эти реакции учитывать 
	// в дополнение к тому, что в угловом смысле ассоциируется с dE/dX.
	mcEndfCrossSectionTable TotalCrossSections;

	// Сечения ядерных реакций в зависимости от энергии падающей частицы
	mcEndfCrossSectionTable NuclearCrossSections;

	std::vector<mcEndfProduct*> Products;

	// Загрузка одного файла сечений
	void Load(const char* fname, const char* ename);

	void Clear();

	void dumpTotalCrossections(std::ostream& os) const;
};
