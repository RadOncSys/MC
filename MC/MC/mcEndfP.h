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

// Class that keeps coss sections for one atomic isotope
class mcEndfP
{
public:
	// Назавание изотопа, включающее атомное имя и атомный вес.
	// Используется как уникальный идентификатор.
	std::string ElementName;

	// Сечения ядерных реакций в зависимости от энергии падающей частицы
	mcEndfCrossSectionTable CrossSections;

	// Загрузка одного файла сечений
	void Load(const char* fname, const char* ename);

	void Clear();
};
