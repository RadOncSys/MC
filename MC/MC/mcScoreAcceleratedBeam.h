// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

/// <summary>
/// Класс анализа пучка ускоренных электронов на уровне мишени.
/// </summary>
class mcScoreAcceleratedBeam : public mcScore
{
public:
	/// <summary>
	/// Конструктор класса анализа пучка ускоренных электронов на уровне мишени.
	/// Радиальная симметрия не подразумевается. Однако, мы пренебрегаем азимутом положения частицы.
	/// Вместо этого используем азимут вектора скорости, так как именно он может приводить
	/// к нарушению радиальной симметри радиационных пучков.
	/// Основанием исключения азимута является малый размер источника и малый вклад в связи с этим в радиационные свойства.
	/// </summary>
	/// <param name="module_name">имя модуля, с которым ассоциирован данный scoring</param>
	/// <param name="nThreads">количество потоков</param>
	/// <param name="ptype">тип частиц, для которых собирается статистика</param>
	/// <param name="ne">количество бинов по энергии</param>
	/// <param name="nr">количество бинов по радиусу</param>
	/// <param name="naxial">количество бинов по азимуту направления частицы</param>
	/// <param name="nthet">количество бинов по углу</param>
	/// <param name="emin">минимальна энергия в МэВ</param>
	/// <param name="emax">максимальная энергия в МэВ</param>
	/// <param name="rmax">максимальный радиус в см</param>
	/// <param name="thetmax">максимальный угол в градусах</param>
	mcScoreAcceleratedBeam(const char* module_name, int nThreads, enum mc_particle_t ptype,
		int ne, int nr, int naxial, int nthet,
		double emin, double emax, double rmax, double thetmax);

	void ScoreFluence(const mcParticle& particle) override;

	void dumpStatistic(ostream&) const override;

	friend ostream& operator << (ostream&, const mcScoreAcceleratedBeam&);

protected:
	enum mc_particle_t ptype_;
	int ne_;
	int nr_;
	int nthet_;
	int naxial_;

	double emin_;
	double emax_;
	double rmax_;
	double thetmax_;

	// Служебные переменные
	double de_;
	double dr_;
	double dthet_;
	double daxial_;

	double skipped_energy_;

	std::vector<double> espec_;				// спектр излучения по всем частицам
	std::vector<double> rspec_;				// распределение энергии по радиусу по всем частицам
	std::vector<double> aspec_;				// распределение энергии по азимуту (вектора скорости) по всем частицам
	std::vector<double> tspec_;				// распределение энергии по направлению по всем частицам

	std::vector<vector<double>> evr_;		// энергетические спектры в зависимости от радиуса
	std::vector<vector<double>> eva_;		// энергетические спектры в зависимости от азимута (вектора скорости)
	std::vector<vector<double>> tva_;		// угловые распределения в зависимости от азимута
	std::vector<vector<double>> wvxy_;		// поток энергии в зависимости отположения
	std::vector<vector<double>> wvvxvy_;	// поток энергии в зависимости от угла
};
