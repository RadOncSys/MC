// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"
#include <vector>

/// <summary>
/// Класс модели фотонного источника с радиальной симметрией, аналогичный прямому сбору частиц (фазовуму пространству).
/// Разнца в представлении частиц (r, Ux, Uy, E), где Ux и Uy - косинусы направления движения в плоскостях
/// прохождения через ось симметрии и текущее положение точки и перпендикулярной ей.
/// Затем, представление в файле целыми числами без знака с масштабированием диапазонов.
/// При таком раскладе на милион частиц потребуется 8 МБ памяти.
/// Учитывая, радиальную симметрию и эффективность размножения с поворотом вокруг оси будем иметь что-то типа 
/// 100 миллионов частиц в транспорте перед коллиматором
/// </summary>
class mcSourceModelRadialDirect : public mcSource
{
public:
	// name		- просто название объекта
	// nThreads - количество обслуживаемых потоков (под каждый свой индекс к общей таблице частиц).
	// z0		- положение источника относительно объекта, к которому он привязан
	mcSourceModelRadialDirect(const char* name, int nThreads, double z0);
	~mcSourceModelRadialDirect();

	void sample(mcParticle& p, mcThread* thread) override;

	// Иницилизация модели (при ее генерации)
	void Init(double maxEnergy, double maxR);

	// Сплитинг устанавливается перед первым самплинго модели.
	// Повторная установка сбросит самплинг в начало.
	void SetSplitting(unsigned nSplits);

	// Наполнение контейнера частиц из массива симуляций.
	void SetParticles(unsigned short* data, unsigned nThreadSize, const std::vector<unsigned>& indexes);
	
	// Доступ к частицам для служебных целей
	const unsigned short* GetParticles() const { return particles_; }
	unsigned GetNoofParticles() const { return noofParticles_; }

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

	void dumpVRML(ostream& os) const override;

	friend std::ostream& operator << (std::ostream&, const mcSourceModelRadialDirect&);

protected:
	// Положение плоскости источника на оси системы
	double z0_;

	// Коэфиициенты, на которые нужно умножить, чтобы получить физические значения параметров
	double energyScale_;
	double rScale_;				// наиболее неопределенная привязка
	double radialAngleScale_;	// компоненты Ux векторов скорости ([-1,1] -> [0,2])
	double azimutAngleScale_;	// компоненты Uy векторов скорости ([-1,1] -> [0,2])

	// Количество частиц в общем контейнере
	unsigned noofParticles_;

	// Количество частиц при расщездении
	unsigned noffsplits_;

	// Индекс текущей частицы при самплинге
	std::vector<unsigned> currentParticleIdx_;

	// Индекс текущей частицы пр сплитинге, который автоматически поддерживается источником.
	// Используем один массив частиц на все потоки, но у каждого потока свои индексы.
	std::vector<unsigned> currentSplitIdx_;

	// Стартовы угол в радианах при сплиттинге частиц
	std::vector<double> splitStartAngle_;
	
	// Зависящий от степени расщепления шаг по углу
	double da_;

	// Коэффициент веса расщепленных частиц
	double dw_;

	// Массив всех частиц в том числе и в режиме симуляции источника.
	// В последнем случае мы разрезаем кусок на равные части для каждого потока.
	// После симуляции склевыем и получаем чистую модель
	unsigned short* particles_;
};
