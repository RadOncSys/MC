// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"
#include "mcHistogramSampler.h"

// Флаг определяющий, будут ли распределения по азимуту апроксимироваться полиномом 
// третьей степени с граничными условиями. 
// В противном случае они будут прямо апроксимироваться гистограммами.
// В процессе длительного изучения оказалось, что при значительном количестве точек 
// эти распределения не подчиняются такому закону. 
// Поэтому, по крайней мере для гамма источников нужно пользоваться прямыми результатами гистограмм.

//#define SRC_USE_ASIMUT_MODEL

/// <summary>
/// Класс модели источника C60 дистанционного аппарата.
/// Данный источник дергает частицы разыгрывая параметры из интегральных диаграм,
/// полученных из реального фазового пространства.
/// </summary>
class mcSourceModelRadialPhsp : public mcSource
{
public:
	mcSourceModelRadialPhsp(const char* name, int nThreads, double z0);

	void sample(mcParticle& p, mcThread* thread) override;
	void sample(mcParticle& p, mcRng& rng) const;

	void dumpVRML(ostream& os) const override;

	friend std::ostream& operator << (std::ostream&,const mcSourceModelRadialPhsp&);

	/// <summary>
	/// Создание образа данных в памяти. Память резервируется внутри функции.
	/// Поэтому вызывающая функция обязана освободить ее с помощью free();
	/// Аргумент size содержит длину выделенной памяти в байтах
	/// </summary>
	void* saveToMemory(int& size);

	/// <summary>
	/// Восстановления модели из образа памяти.
	/// </summary>
	void readFromMemory(void* buffer);

	/// <summary>
	/// Создание самплеров из матрицы симуляции.
	/// </summary>
	void createFromDistribution(int isxray, int ne, int nr, int nt, int na,
		double emax, double r2max, double tmax, double focus,
		const double* data);

	static double sampleAzimut(const vector<double>& a, double x);
	static void fitAsimutDistribution(int np, const double* data, double& a, double& b);

protected:
	/// <summary>
	/// Получение положения и углов для отдельно взятого эннергетического бина.
	/// </summary>
	void sampleForEnergy(int eidx, double& r, double& t, double& a, double rnd0, double rnd1, double rnd2) const;

	/// <summary>
	/// Создание массивов самплеров.
	/// </summary>
	void createHistograms(int ne, int nr, int nt, int na);

protected:
	int isxray_;

	// Положение плоскости источника на оси системы
	double z0_;

	int ne_;
	int nr_;
	int nthet_;
	int naxial_;

	double focus_;
	double emax_;
	double rmax_;
	double thetmax_abs_;
	double thetmax_;

	// Служебные геометрические переменные
	double de_;
	double dr_;
	double dthet_abs_;
	double dthet_;
	double daxial_;

	// Веса линий фотонов
	double we117_;
	double we133_;
	double wexscale_;

	// Гистограммы для самплинга положения и углов
	mcHistogramSampler xehistogram_;
	vector<mcHistogramSampler> rhistograms_;
	vector<vector<mcHistogramSampler>> thistograms_;

#ifdef SRC_USE_ASIMUT_MODEL
	vector<vector<vector<vector<double>>>> apars_;
#else
	vector<vector<vector<mcHistogramSampler>>> ahistograms_;
#endif
};
