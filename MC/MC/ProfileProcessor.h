// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <vector>
#include <memory>

// Класс анализа дозовых профилей.
class ProfileProcessor
{
public:
	/// <summary>
	/// Определение параметров поперечного профиля.
	/// </summary>
	/// <param name="x">координаты точек профиля</param>
	/// <param name="d">значения доз</param>
	/// <param name="width">расчитанная ширина профиля по полувысоте</param>
	/// <param name="penumbra">расчитанный размер полутении</param>
	/// <param name="penumbra">дозиметрический индекс, соответствующий интегралу провиля от 80% плюс 3 см</param>
	static void ProfileParameters(const std::vector<double>& x, const std::vector<double>& d,
		std::vector<double>& dd, double& width, double& penumbra, double& dinde);

	/// <summary>
	/// Сглаживание дозовых распределений круглых полей в FanLine RZ геометрии.
	/// Идея сглаживания в по радиусу в окрестности центральной оси апроксимацией параболой.
	/// Учитываются точки на расстоянии не более половины размера поля и не ближе к его границе чем 1.5 см.
	/// При этом в весе точек присутствует статистический коэффициент погрешности пропорциональный объему детектора.
	/// По глубине распределения за исключением области build-up меняются медленно.
	/// Поэтому сглаживаем их одномерным Savitzky–Golay высокой степени.
	/// </summary>
	/// <param name="srcMatrix">исходная шумная матрица точек (данные упакованы по радиусу глубина за глубиной)</param>
	/// <param name="nr">размер матрицы по радиусу</param>
	/// <param name="nz">размер матрицы по глубине</param>
	/// <param name="rstep">шаг по радиусу в см (для определения пределов сглаживания по радиусу)</param>
	static std::shared_ptr<std::vector<std::vector<double>>>
		SmoothFanRZ(const std::vector<std::vector<double>>& srcMatrix, unsigned nr, unsigned nz, double rstep);

	// Сглаживание начального участка радиально симметричного профиля параболой
	// nr - количествао точек от начала, подверженных сглаживанию.
	static void SmoothRZProfile(std::vector<double>&, unsigned nr);

	// Одномерное сглаживание Savitzky–Golay с различными степенями.
	static void SmoothSG1D(std::vector<double>& p, unsigned m, unsigned nl, unsigned nr);

	// Двумерное сглаживание Savitzky–Golay.
	static std::unique_ptr<std::vector<std::vector<double>>> SmoothSG2D(const std::vector<std::vector<double>>& data);
};
