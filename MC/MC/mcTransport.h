// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcParticle.h"
#include "mcObj.h"
#include "mcScore.h"
#include <vector>

class mcMedia;
class mcRng;
class mcScore;
class mcMediumXE;

#ifndef NNEG
#define  NNEG(a)    ((a) >= 0  ? (a) : -(a))
#endif

#ifndef ZORP
#define  ZORP(a)    ((a) >= 0  ? (a) : 0)
#endif

#ifndef MAX
#define  MAX(a, b)  ((a) > (b) ? (a) :  (b))
#endif

#ifndef MIN
#define  MIN(a, b)  ((a) < (b) ? (a) :  (b))
#endif

#define MC_EPSILON     1e-6 // используется в борьбе с ошибками округления

enum mc_move_result_t { MCMR_INTERUCT = 0, MCMR_EXIT, MCMR_CONTINUE, MCMR_DISCARGE };

/// <summary>
/// Базовый класс транспорта, являющийся ключевым в дизайне
/// всего модуля симуляции методом Монте-Карло.
///
/// Система транспорта предполагает, что в мировой системе координат все объекты
/// выстроены вдоль оси Z, причем, координаты Z этих объектов не пересекаются.
/// Это позволяет принципиально упростить решение двух задач.
///
/// Во-первых можно выстроить линейную цепочку объектов / транспортов.
/// При покидании текущего объекта транспорт передается в предыдущий или следующий 
/// объект в зависимости от направления движения вдоль оси Z.
///
/// Во-вторых, система координат частицы преобразуется в систему координат объекта,
/// что упрощает расчеты геометрии и scoring.
/// 
/// GG 20150227
/// Базовый класс модулей транспорта объектов, вложенных друг в друга как Русская матрешка.
/// Класс поддерживает интерфейс транспорта и дополнительно имеет ссылки на вложенный и окружающий объект.
/// Расстояние до границы определяется как расстояние до внешней поверхности и до вложенного объекта если таковой имеется.
/// Таким образом, скажем, свферическое кольцо можно рассматривать как сферу, внутри которой находится другая сфера.
/// Решение о переходе в другой объект принимается на основании того, какя поверхность будет пересечена.
/// В описании сцены группа вложенных объектов помечается именно как группа со своими глобальными координатами.
/// Система же координат каждого объекта указывается в XML файле относительно группы.
/// Однако, во избежание лишнего посредника объекты группы получают систему координат сцены путем объединения с преобразованиями группы.
/// Базовый класс траспорта предоставляет функцию определения объекта для последующего транспорта подобно своему базовумум классу.
/// Разница в том, что дополнительно анализируются вложенные и окружающие объекты.
/// </summary>
class mcTransport : public mcObj
{
public:
	mcTransport();
	mcTransport(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransport(void);

	//
	// Инициализация
	//

	/// <summary>
	/// Позиционирование объекта в мировой системе координат.
	/// Это единственное место, где устанавливаются матрицы преобразования координат.
	/// </summary>
	void setPosition(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);

	/// <summary>
	/// Ссылка на объект, содержащий все характеристики сред и ссылки на физические расчеты
	/// </summary>
	void setMediaRef(const mcMedia* media);

	/// <summary>
	/// Функция, предназначенная для установки среды транспорта, состоящего из одного монолитного тела.
	/// Стандартная функция перемещения частицы использует указанную среду как среду объекта.
	/// Вне среды тела предполагается вакуум. 
	/// Если потребуется учитывать воздух, между объектами можно размещать слои воздуха.
	/// </summary>
	void setMediumMono(const char* mname);

	/// <summary>
	/// Аналогично плотность этого объекта, по сравнению с физической плотностью в нормальных условиях
	/// </summary>
	void setMediumMono(double d) { defdensity_ = d; }

	/// <summary>
	/// Коррекция систем преобразования координат, изначально созданная для поддержки 
	/// задания положения объектов относительно группы
	/// </summary>
	void MoveToCoordinateSystem(const geomMatrix3D& m);

	//
	// Транспорт
	//

	/// <summary>
	/// Транспорт частицы прилетающей извне объекта данного транспорта начинается 
	/// с функции beginTransport, в которой сначала происходит преобразование координат
	/// из мировой системы в локальную и запускается физика.
	/// Транспорт создает копию частицы, поэтому ее можно объявить константой.
	/// </summary>
	virtual void beginTransport(mcParticle& p);

	/// <summary>
	/// Транспорт частицы рождающейся внутри объекта данного транспорта начинается. 
	/// В отличие от функции beginTransport частице сразу назначается регион и она не переносится на поверхность.
	/// </summary>
	virtual void beginTransportInside(mcParticle& p);

	/// <summary>
	/// Когда частица покидает объект, ее координаты и направление преобразуются обратно в мировую систему.
	/// Затем, становится ясно, какому транспорту передать частицу.
	/// Речь идет о текущей активной частице.
	/// </summary>
	virtual void endTransport(mcParticle* particle);

	/// <summary>
	// Симуляция текущей частицы в указанном потоке. 
	/// Симуляция не замещается в наследующих классах, так как алгоритм является универсальным.
	/// Индивидуализация осуществляется через виртуальные методы.
	/// </summary>
	static void simulate(mcThread* thread);

	/// <summary>
	/// Перемещаем частицы на заданное количество (хранящееся в самой частице) средних пробегов.
	/// Конкретный класс определяет свойства регионов и отслеживает ссылки.
	/// Стандартная функция использует чисто геометрические методы, предполагая монолитность
	/// объекта (т.е., состоящего из одного региона) и отсутствие скоринга.
	/// В аргументах возвращается реальный шаг и переданная на пути этергия
	/// (непрерывные потери заряженных частиц).
	/// </summary>
	virtual mc_move_result_t moveParticle(mcParticle* particle, double& step, double& edep);

	//
	// Scoring
	//

	/// <summary>
	/// Класс транспорта уничтожает объект score в деструкторе.
	/// Это важно помнить и в запрашивающей программе score 
	/// следует создавать по new и не вызывать delete по окончании!!!
	/// </summary>
	void setScore(mcScore* score);
	mcScore* getScore() { return score_; }

	/// <summary>
	/// Устанавливает скоринг, который только накапливает энергию, выделившуюся в объете транспорта.
	/// </summary>
	void setDefaultScore(int nThreads);

	/// <summary>
	/// Добавление объекта (региона) в группу.
	/// Объектом является любой объект транспорта.
	/// Объект должен быть сосздан с системой координат относительно смистемы координат объекта группы.
	/// Функция добавления сама произведет необходимые преобразования по переведению системы координат относительно мировой.
	/// </summary>
	void addRegion(mcTransport* region);

	//
	// Вспомогательные методы
	//
	void setPreviousTransport(mcTransport* t) { previousTransport_ = t; }
	void setNextTransport(mcTransport* t) { nextTransport_ = t; if (t != nullptr) t->setPreviousTransport(this); }
	const mcTransport* getPreviousTransport() const { return previousTransport_; }
	const mcTransport* getNextTransport() const { return nextTransport_; }

	void setInternalTransport(mcTransport* t);
	void setExternalTransport(mcTransport* t);
	mcTransport* getInternalTransport() const { return internalTransport_; }
	mcTransport* getExternalTransport() const { return externalTransport_; }

	// Запрос вложенного транспорта по имени
	virtual mcTransport* getInternalTransportByName(const char* name);

	void setStamp(int stamp) { stamp_ = stamp; }

	// Преобразования координат
	const geomMatrix3D& MW2T() { return mwtot_; }
	const geomMatrix3D& MT2W() { return mttow_; }

	// Разыгрывание пути
	static double HowManyMFPs(mcRng& rng);

	double etotal() const;

	static void implementException();

	virtual void dump(ostream& os) const;
	virtual void dumpVRML(ostream& os) const;

	void dumpVRMLRing(ostream& os, double r1, double r2, double z, bool normPositive, double x0 = 0, double y0 = 0) const;
	void dumpVRMLSemiCircle(ostream& os, double r, double z, bool normPositive, const geomMatrix3D& M) const;
	void dumpVRMLCylinderSide(ostream& os, double r, double z1, double z2, bool normPositive, double x0 = 0, double y0 = 0) const;
	void dumpVRMLCylinderSemiSide(ostream& os, double r, double h, const geomMatrix3D& M) const;
	void dumpVRMLConicalCylinderSide(ostream& os, double r, double z1, double z2, double f, bool normPositive, double x0 = 0, double y0 = 0) const;
	void dumpVRMLCylinderRing(ostream& os, double r1, double r2, double z1, double z2) const;
	void dumpVRMLConicalRing(ostream& os, double r1, double r2, double z1, double z2, double f) const;
	void dumpVRMLCylinderWithConicalHole(ostream& os, double r1, double r2, double z1, double z2, double f) const;
	void dumpVRMLRectangleRing(ostream& os, double x1, double x2, double y1, double y2, double d, double h) const;
	void dumpVRMLCylinder(ostream& os, double r, double z1, double z2, double x0 = 0, double y0 = 0) const;
	void dumpVRMLPrism(ostream& os, double ax, double ay, double az) const;
	void dumpVRMLPolygonCircle(ostream& os, const std::vector<double>& pz, const std::vector<double>& pr) const;

	//
	// Геометрия (методы нужны для поддержки стандартной функции moveParticle)
	//
	virtual double getDistanceInside(mcParticle& p) const;
	virtual double getDistanceOutside(mcParticle& p) const;
	virtual double getDNearInside(const geomVector3D& p) const;

	//
	// Auxilary
	//
	short getDefMedIdx() const { return defmedidx_; }
	double getDefDensity() const { return defdensity_; }

	// GG 20140927 Механизм регулирования ускорения транспорта через пороги уничтожения частиц на уровне объекта транспорта.
	// Частица имеет ссылку на объет транспорта и это можно использовать не меняя инфраструктуру.
	// По умолчанию пороги устанавливаются в 0.
	// В этом случае используется порог из Media как и раньше.
	// Внутри транспорта пороги электронов используются по одному разу в физике электронов и позитронов и три раза в физике фотонов.
	// В физике фотонов это единственное место, связанное с регулирование фото электронов и трогать его не надо.
	// В зараяженных частицах это только функция Discharge.
	double transCutoff_phot; // Energy cutoff for photon transport
	double transCutoff_elec; // Energy cutoff for electron transport

protected:
	virtual mcMediumXE* getParticleMedium();
	virtual double getRegionRelDensity();

protected:
	mcScore* score_;
	mcTransport* previousTransport_;
	mcTransport* nextTransport_;
	mcTransport* internalTransport_;
	mcTransport* externalTransport_;

	geomMatrix3D mwtot_;	// преобразование из мировой системы в систему объекта
	geomMatrix3D mttow_;	// преобразование из системы объекта в мировую систему
	geomMatrix3D mtoe_;		// преобразование из данного объекта во вложенный

	const mcMedia* media_;
	short defmedidx_;
	double defdensity_;

	// Индивидуальная печать данного объекта транспорта
	int stamp_;

	bool isMultiRegions_;
	std::vector<mcTransport*> regions_;
};
