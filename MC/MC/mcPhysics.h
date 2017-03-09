// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

class mcParticle;
class mcMedium;
class geomVector3D;

// Базовый класс моделирования физических процессов
class mcPhysics
{
public:
	mcPhysics(void);
	virtual ~mcPhysics(void);

	// Расчет средней длины пробега в среде.
	virtual double MeanFreePath(double ke, const mcMedium& med, double dens) const = 0;

	// Разыгрывание акта взаимодействия.
	// Подразумевается, что указатель на часицу является указателем на элемент массива (стэка).
	// Если в результате взаимодействия вместо одной частицы образуется несколько, лишние
	// поместятся в стэк на основании указанного адреса текущего курсора.
	// Возвращается выделившаяся в точке в результате взаимодействия энергия.
	virtual double DoInterruction(mcParticle* p, const mcMedium* med) const = 0;

	// Перемещение частицы на заданное расстояние. Возвращает переданную на пути энергию.
	// Для фотонов это просто изменение координаты частицы.
	// В переменной step возвращается реальный шаг. Для заряженных частиц он будет отличаться от заданного,
	// если шаг моделирования непрерывных потерь окажется меньше.
	// В последнем случае транспорт определяет, что дискретное событие разыгрывать не следует.
	virtual double TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const = 0;

	// Проверка частицы на предмет энергии ниже критической и при необходимости ее уничтожение.
	// Выделяющаяся энергия помещается в edep.
	// Уничтожение именно здесь, так как возможны процессы типа аннигиляции, что известно только физике.
	virtual bool Discarge(mcParticle* p, const mcMedium& med, double& edep) const = 0;

	// Physics utilities
	static void GetRandomPhi(double rnum, double* cosPhi, double* sinPhi);
	static void ChangeDirection(double cosTheta, double sinTheta, double cosPhi, double sinPhi, geomVector3D& u);
	static void GoInRandomDirection(double rnum1, double rnum2, geomVector3D& u);

	static mcParticle* DuplicateParticle(mcParticle* p);
	static void DiscardParticle(mcParticle* p);
};
