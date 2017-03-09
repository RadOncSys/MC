#include ".\mcgeomrectside.h"
#include <float.h>

mcGeomRectSide::
mcGeomRectSide(const geomVector3D& p
	, const geomVector3D& Vx
	, const geomVector3D& Vy
	, double ax, double ay)
	:P_(p)
	, Vx_(Vx)
	, Vy_(Vy)
	, ax_(ax)
	, ay_(ay)
{
	Vx_.normalize();
	Vy_.normalize();
	N_ = Vx_^Vy_;
}

double mcGeomRectSide::getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const
{
	// Высота точки над плоскостью, в которой определен прямоугольник
	double h = (P_ - p)*N_;

	// Наклон вектора направления к плоскости
	double cosn = v*N_;
	if (fabs(cosn) < DBL_EPSILON)
		return DBL_MAX; // параллельно плоскости

	if (fabs(h) < DBL_EPSILON && ((inside && cosn < 0) || (!inside && cosn > 0)))
		return DBL_MAX; // движение от поверхности

	// Точка пересечения вектора направления с плоскостью
	double distance = h / cosn;
	if (distance < 0)
		return DBL_MAX; // движение от плоскости
	geomVector3D pc = p + (v*distance);

	// Координаты точки проекции в системе прямоугольника
	geomVector3D pr = pc - P_;
	double x = pr*Vx_;
	double y = pr*Vy_;

	// Если точка пересечения не в рамке, то расстояние бесконечное
	if (x >= 0 && x <= ax_ && y >= 0 && y <= ay_)
		return distance;
	else
		return DBL_MAX; // в нужном направлении но мимо прямоугольника
}

double mcGeomRectSide::getDNear(const geomVector3D& p) const
{
	// Высота точки над плоскостью, в которой определен прямоугольник
	double h = (p - P_)*N_;

	// Проекция точки на плоскость прямоугольника (по нормали)
	geomVector3D pp = p + (N_*h);

	// Координаты точки проекции в системе прямоугольника
	geomVector3D pr = pp - P_;
	double x = pr*Vx_;
	double y = pr*Vy_;

	// Если проекция в рамке прямоугольника, то минимальным расстоянием будет высота
	if (x >= 0 && x <= ax_ && y >= 0 && y <= ay_)
		return fabs(h);

	// Перебор областей
	if (x < 0)
	{
		if (y < 0) return sqrt(x*x + y*y + h*h);
		else if (y < ay_) return sqrt(x*x + h*h);
		else return sqrt(x*x + (y - ay_)*(y - ay_) + h*h);
	}
	else if (x > ax_)
	{
		if (y < 0) return sqrt((x - ax_)*(x - ax_) + y*y + h*h);
		else if (y < ay_) return sqrt((x - ax_)*(x - ax_) + h*h);
		else return sqrt((x - ax_)*(x - ax_) + (y - ay_)*(y - ay_) + h*h);
	}
	else
	{
		if (y < 0) return sqrt(y*y + h*h);
		else return sqrt((y - ay_)*(y - ay_) + h*h);
	}
}
