﻿// Radiation Oncology Monte Carlo open source project
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
#include "mcRng.h"

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

	// Парсинг целого числа
	static int iStrCrop(const char* s, int n) {
		std::string s1 = s;
		s1.erase(n);
		int f = std::stoi(s1);
		return f;
	}
};

// Crossections per istope per incident particle energy
class mcEndfCrossSectionTable
{
public:
	mcEndfCrossSectionTable() {
		isEmpty = true;
		MT = -1;
	}

	void Load(std::istream& is);
	void dump(std::ostream& os) const;

	double get_lambda(double kE, double rho, double A);

	double get_sigma(double kE) const;

	//Статус, показывающий не пуста ли таблица
	bool isEmpty;

	// Количество пар энергия падающей частицы / сечение
	std::vector<int> npoints;

	// Количество видов интерполяции
	// Временно предполагаем, что мы не столкнемся со 
	// множесвенностью интерполяций в пределах одной таблицы.
	// TODO: Проверить эту гипотезу (в коде стоит exception на этот случай)
	int ninterpolations;
	 
	// Тип интерполяции.
	// TODO: если обнаружится потребность в поддержке 
	// множества типов интерполяций переделать в массив
	std::vector<int> interpolationType;

	short MT;
	
	// Точки
	std::vector<double> Energies;
	std::vector<double> Values;
};

// Класс для чтения MF=6 MT=5 

class mcEndfEANuclearCrossSectionTable
{
public:
	mcEndfEANuclearCrossSectionTable() {
		isEmpty = true;
	}
	void mLoad(std::istream& is);
	
	void Load(std::istream& is, int LAW);

	void dump(std::ostream& os) const;

	//Розыгрыш мультиплетности с интерполяцией
	int playMulti(double kE, mcRng& rng) const;

	//Розыгрыш энергии вылетающей частицы
	double playE(double kE, int &keIN, int &eoutID, mcRng& rng) const;

	//Розыгрыш f_0 и r
	double** playpar(mcRng& rng, double kE, int LAW);

	//Розыгрыш косинуса угла рассеяния
	double playmu(double kE, int LAW, int keIN, int eoutID, int ptype, mcRng& rng) const;

	//Интерполяция f_0 для пары энергия-энергия вылета
	double getf_0(int IN, double Eout);

	double integrate_f0(mcRng& rng, double kE);

	int ZA_nucl;

	double AWR_nucl;

	bool isEmpty;

	//Интерполяция 

	// Количество пар энергия падающей частицы / мультиплетность
	int n_energypoints;

	// Количество энергий вылета проодуктов в распределении (NEP в ENDF [Chapter 6])
	int npoints_out;

	//LANG (= 1 - представление Лежандра, = 2 - представление Кальбаха-Манна)
	std::vector<int> LANG;

	int iLang;

	//Количество угловых параметров
	int NA;

	// Количество видов интерполяции
	// Временно предполагаем, что мы не столкнемся со 
	// множесвенностью интерполяций в пределах одной таблицы.
	// TODO: Проверить эту гипотезу (в коде стоит exception на этот случай)
	int ninterpolations;

	// Тип интерполяции.
	// TODO: если обнаружится потребность в поддержке 
	// множества типов интерполяций переделать в массив
	int interpolationType;
	
	//Трехмерный вектор с энерго-угловыми параметрами
	std::vector<std::vector<std::vector<double>>> EA_par;

	std::vector<double> EA_Epoints;
	
	//LAW = 1:
	//EA_par parameters for Kalbach-Mann:
	//EA_par[i][j][k], where i - identify incident energy
	//						 j - identify outer energy
	//					and  k - identify corresponding parameter
	//Exactly if NA = 1 then [i][j][0] keeps E_out (from incedent E_i to E_out)
	//						 [i][j][1] keeps f_0 (total emission probability from E_i to E_out)
	//						 [i][j][2] keeps "r" (special parameter)
	// if NA = 2 then added  [i][j][3] keeps "a" (second special parameter)
	//For Legandre representation i, j and k have the same meaning and
	//						 [i][j][0] keeps E_out
	//						 [i][j][1] keeps f_0
	//						 [i][j][2] keeps f_1
	//							...
	//						 [i][j][NA+1] keeps f_NA
	//LAW = 2:
	// LANG = 12:
	// EA_par[i][j][0] keeps incident energy of i-th proton (incedent particle)
	// EA_par[i][j][1] keeps cosine of scattering j-th angle
	// EA_par[i][j][2] keeps p(mu) - differential probability to scatter at this j-th angle (for LANG = 12)
	// LANG = 0: (Legandre coefficients)
	// EA_par[i][0][0] keeps incident energy of i-th proton (incedent particle)
	// EA_par[i][1][0] keeps a_1
	// EA_par[i][2][0] keeps a_2
	//			...
	// EA_par[i][NL][0] keeps a_NL
	// 
	//

	// Точки
	std::vector<double> Energies;
	std::vector<double> Multiplicities;
};

enum particle_type { neutron = 0, proton, deutron, triton, alpha, recoils, gammas, electron };

class mcEndfProduct
{
public:

	mcEndfProduct();

	~mcEndfProduct();

	//Type of product
	particle_type product_type;

	int ZAP;

	double AWP;

	//Закон представления распределения
	int LAW;

	// Энерго-угловые сечения в зависимости от энергии налетающих протонов
	std::vector<std::shared_ptr<mcEndfEANuclearCrossSectionTable>> EANuclearCrossSections;
	
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

	// Сечения реакций (p,n) MF = 3 MT = 50
	mcEndfCrossSectionTable Neutron0CrossSection;

	// Сечения реакций (p,n) MF = 3 MT = 51
	mcEndfCrossSectionTable Neutron1CrossSection;

	// Сечения реакций (p,n) MF = 3 MT = 52
	mcEndfCrossSectionTable Neutron2CrossSection;					//СДЕЛАТЬ VECTOR ДЛЯ MT = 50 - 90

	// Сечения реакций (p,n) MF = 3 MT = 53
	mcEndfCrossSectionTable Neutron3CrossSection;

	// Сечения реакций (p,n) MF = 3 MT = 54
	mcEndfCrossSectionTable Neutron4CrossSection;

	// Сечения реакций (p,n) MF = 3 MT = 55
	mcEndfCrossSectionTable Neutron5CrossSection;

	std::vector<mcEndfProduct*> Products;

	//MT = 50, 51...; MF = 6;
	std::vector<mcEndfProduct*> EmittedNeutrons;

	// Загрузка одного файла сечений
	void Load(const char* fname, const char* ename);

	void Clear();

	void dumpTotalCrossections(std::ostream& os) const;
};

// Класс для чтения MF=4

class mcEndfAngular
{
public:
	mcEndfAngular() {
		isEmpty = true;
		NE1 = 0;
		NE2 = 0;
		MT = -1;
	}

	void Load(std::istream& is);

	short MT;

	bool isEmpty;

	int ZA;

	double AWR;

	// Вид представления угловых распределений в секции
	short LTT;

	// 1 - все распределения изотропны, 0 - нет
	short LI;

	// Система отсчета, в которой представлены распределения
	short LCT;

	int NE1, NE2;
	
	//Распределения по Лежандру
	std::vector<double> LEnergies;
	std::vector<std::vector<double>> LValues;

	//Табличные распределения
	std::vector<double>TEnergies;
	std::vector<std::vector<double>> Cosines;
	std::vector<std::vector<double>> TValues;

};

class mcEndfN
{
public:
	mcEndfN();

	~mcEndfN();

	// Назавание изотопа, включающее атомное имя и атомный вес.
	// Используется как уникальный идентификатор.
	std::string ElementName;

	mcEndfCrossSectionTable TotalCrossSections;

	mcEndfCrossSectionTable ElasticCrossSections;

	mcEndfCrossSectionTable InelasticCrossSections;

	// Сечения ядерных реакций в зависимости от энергии падающей частицы
	mcEndfCrossSectionTable NuclearCrossSections;

	// MF = 4, MT = 2
	mcEndfAngular nElasticAngular;
	
	std::vector<mcEndfAngular*> inelasticLevelsAng;

	// Сечения реакций (p,n) MF = 3 MT = 51-90;
	std::vector<mcEndfCrossSectionTable*> nInelasticCS;

	std::vector<mcEndfProduct*> Products;

	//MT = 91; MF = 6;
	std::vector<mcEndfProduct*> nInelasticContin;

	// Загрузка одного файла сечений
	void Load(const char* fname, const char* ename);
};