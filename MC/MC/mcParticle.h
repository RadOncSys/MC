// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcRegionReference.h"
#include "../geometry/vec3d.h"

class mcScoreTrack;
class mcThread;
class mcTransport;

// Типы частиц введены в качестве индексов объектов транспорта
enum mc_particle_t { MCP_PHOTON = 0, MCP_NEGATRON, MCP_POSITRON, MCP_PROTON, MCP_NEUTRON, MCP_NTYPES };

class mcParticle
{
public:
	mcParticle(void);
	mcParticle(mc_particle_t pt, int pq, double pke, const geomVector3D& pp, const geomVector3D& pu);
	mcParticle(const mcParticle& p);
	~mcParticle(void);

public:
	enum mc_particle_t t;
	int q;
	double ke;
	geomVector3D p;
	geomVector3D plast;	// Точка последнего взаимодействия
	geomVector3D u;
	double dnear;
	mcRegionReference region;
	double weight;

	// Расстояние до следующего взаимодействия в единицах длин среднего пробега
	double mfps;

	// Относительная плотность региона, в котором находится частица в данный момент
	double regDensityRatio;

	// Указатель на транспортный объект,
	// в котором частица находится в данный момент
	mcTransport* transport_;

	// Указатель на транспортный объект с минимальным расстоянием в случае возможности 
	// передачи управления модулю группы
	mcTransport* transportNearest_;

	// Место рождения частицы
	//mcRegionReference reg_born;

	// Место последнего взаимодействия
	//mcRegionReference reg_last_scatter;

	// Указатель на объект, собирающий перемещения частиц между модулями
	mcScoreTrack* trackScore_;

	// Указатель на объект, содержащий специфичные для потока данные.
	// Поскольку именно над частицей производятся манипуляции,
	// нуждающиеся в использовании специфичных для потока данных,
	// то именно через нее посредством объекта потока передается вся информация, 
	// включая генераторы случайных чисел и идентификатор потока.
	mcThread* thread_;

	// Объекты транспорта могут ставить печать в момент 
	// рождения частицы. Рассеянная частица означает ее рождение
	int regionBirth;

	// Объекты транспорта могут ставить печать в момент 
	// BeginTransport обозначая, что частица в них находилась.
	int regionFlags;

	// Тип пересекаемой поверхности (внешняя или внутрення) в embedded transport.
	// Не определенное значение полезно для контроля за тем, что соответствующие 
	// функции вычисления расстояний учитывают необходимость установки данного флага.
	enum temb_shit_t : short { Undefined = 0, External, Internal };
	temb_shit_t exitSurface_;
};
