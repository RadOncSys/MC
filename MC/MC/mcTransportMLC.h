// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"
#include <vector>

// Класс MLC.
// Для начала без лепестков.
// Только прямоугольные поля правдоподобной формы
class mcTransportMLC : public mcTransport
{
public:
	mcTransportMLC(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		double focus, double r, double h, double w, double l, 
		double fsx1, double fsx2, double fsy1, double fsy2);
	virtual ~mcTransportMLC(void);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	// Перед обращением к расстояниям до блока координаты и направления частицы 
	// приводятся к системе блока, заключающейся в контакте скругленной части с линией X=0
	// при расположении самого блока слева.
	double getDistanceToBlockInside(const geomVector3D& p, const geomVector3D& u, double y1, double y2, double xx) const;

	// Вычисление расстояния до MLC при нахождении за внешним периметром
	double getDistanceToMlcOutsideField(const geomVector3D& p, const geomVector3D& u) const;

	// Вычисление расстояния до MLC при нахождении внутри радиационного поля
	double getDistanceToMlcInsideField(const geomVector3D& p, const geomVector3D& u) const;

	// Вычисление расстояния до MLC при нахождении внутри отдельно взятой части радиационного поля.
	// Здесь мы возвращаем не бесконечность, если частица попадает в торец лепестков блока или
	// пересекает боковую поверхность (XZ). Какой из последних случаев должна определять вызывающая функция.
	double getDistanceToMlcInsideFieldBlock(const geomVector3D& p, const geomVector3D& u, double x1, double x2, double y1, double y2) const;

	// Вычисление расстояния до конуса, образованного снаружи МЛК утопленной центральной частью.
	// Если пересечение с боковыи поверхностями за пределами слоя МЛК, то возвращается бесконечность.
	double getDistanceToMlcRectangleConeInside(const geomVector3D& p, const geomVector3D& u, double x1, double x2, double y1, double y2) const;

	// Проверка лежит ли точка на поверхности коллиматора в плоскости Z = const
	bool isInZPlaneHit(double x, double y) const;

	// Предполагается, что частица находится внутри слоя коллиматора.
	// Определяется, находится ли она за пределами внешнего периметра.
	bool isOutsideMLC(double x, double y) const;

	// Предполагается, что точка на одной из граней XZ и в слое MLC (по Z).
	// Определяем, лежит ли точка на лепестке.
	bool isOnOutsideLeaf(const geomVector3D& p) const;

	// Определяем, лежит ли точка на лепестке.
	// Только левый лепесток, так как предполагаетс, что все уже развернуто в левую сторону.
	bool isOnOnTheLeftLeaf(const geomVector3D& p, double y0, double xx) const;

	void dumpVRMLMLC(ostream& os) const;

	double focus_;	// расстояние от источника до дальней поверхност коллиматора
	double r_;		// радиус кривизны торца лепестков
	double h_;		// высота лепестка
	double w_;		// ширина банка лепестков
	double l_;		// длина отдельного лепестка
	double fsx1_;	// размер поля
	double fsx2_;	// (физические координаты наиболее выступающей части лепестка,
	double fsy1_;	// т.е. без коррекции для обеспечения нужного размера поля
	double fsy2_;	// на уровне изоцентра)

	// Вспомогательные переменные для ускорения
	double r2_;
	double dx_;		// отклонение скругленной части от наиболее выступающей по X
	double sfy_;	// масштабный коэффициент перевода по Y с передней на заднюю поверхность

	geomVector3D pymin_, pymax_;	// служебные вектора (0, +/-w/2, 0)
	geomVector3D py1_, py2_;		// служебные вектора (0, +fsy1/fsy2, 0)
	geomVector3D zxn1_, zxn2_;		// нормали к граням XZ
};
