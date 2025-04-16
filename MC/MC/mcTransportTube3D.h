// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// TODO: Класс стартовал с существующего транспорта mcETransportConvexPolygonCircle.
// Отличие только в том, что ось типа радиально симметричной фигуры изгибается.
// Далее нужно заменить функционал mcETransportConvexPolygonCircle 
// на правильный для mcTransportTube3D

// Класс данных отдельного сегмента
struct TubeSection
{
	geomVector3D P0;	// Положение центра в системе трубки
	geomVector3D V0;	// Направление оси сегмента
	geomVector3D N0;	// Нормаль к стартовой секущей сегмент плоскости
	geomVector3D NT;	// Нормаль к стартовой секущей сегмент плоскости в системе трубки
	geomVector3D NT1;	// Нормаль к плоскости второго торца в системе трубки
	double R;			// Внешний радиус
	double D;			// Расстояние между точками P0 данного и следующего сегментов

	// Преобразование координат из системы эндостата 
	// в стартовую пллоскость сечения сегмента
	geomMatrix3D ME2SPlane;

	// Преобразование координат из системы эндостата в цилиндр сегмента.
	geomMatrix3D ME2STube;
};

// Класс транспорта в изгибающейся полнотелой трубке типа эндостата брахитерапии.
// При симуляции реального эндостато нужно создавать два объекта этого типа.
// Один (воздух) вложенный в другой (материал эндостата).
class mcTransportTube3D : public mcTransport
{
public:
	mcTransportTube3D(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, 
		const std::vector<geomVector3D>& p, const std::vector<double>& r);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	// Нахождение сегмента в котором находится точка p.
	// Возможны значения -1 и nsegments, 
	// если точка до стартового или после последнего сегмента.
	int findSegment(const geomVector3D& p) const;

	// Расстояние до боковой стенки сегмента трубки.
	// Предполагается, что до вызова уже проведена проверка столкновения с торцевыми плоскостями 
	// и установлено, что пересечение с боковой стенкой сегмента возможно.
	// Боковая стенка не цилиндр и не конус из-за различий радиусов торцов
	// и наклона плоскости вторго торца.
	double segmentTubeDistanceInside(int idx, const geomVector3D& p, const geomVector3D& u) const;
	double segmentTubeDistanceOutside(int idx, const geomVector3D& p, const geomVector3D& u) const;

	// Расстояние до точки пересечения поверхности, соответствующей точке пересечения траекторией (pc)
	double getSurfaceR(int idx, const geomVector3D& pc) const;
	double getSurfaceRCut(int idx, const geomVector3D& pc) const;

	std::unique_ptr<std::vector<TubeSection>> segments_;
};
