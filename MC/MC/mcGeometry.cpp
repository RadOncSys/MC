#include "mcGeometry.h"
#include "mcDefs.h"
#include "../geometry/vec3d.h"
#include <float.h>

double mcGeometry::getDistanceToInfiniteCylinderInside(const geomVector3D& p, const geomVector3D& v, double r)
{
	double cd = 0;
	double a = v.sqLengthXY();
	if (a < DBL_EPSILON) return DBL_MAX;		// Параллельно оси
	double b = (p.x()*v.x() + p.y()*v.y()) / a;
	double c = (p.sqLengthXY() - r*r) / a;
	double det = b*b - c;
	if (det > 0)
		cd = -b + sqrt(det);  // when points is inside, only one positive solution exists
	return (cd < 0 ? 0 : cd);
}

double mcGeometry::getDistanceToInfiniteCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r)
{
	double b = (p.x()*v.x() + p.y()*v.y());
	if (b >= 0) return DBL_MAX;
	double a = v.sqLengthXY();
	if (a < DBL_EPSILON) a = DBL_EPSILON;
	b /= a;
	double c = (p.sqLengthXY() - r*r) / a;
	double det = b*b - c;
	if (det <= 0) return DBL_MAX;
	double cd = (-b - sqrt(det));
	if (cd < 0) return DBL_MAX;
	else return cd;
}

double mcGeometry::getDistanceToCylinderInside(const geomVector3D& p, const geomVector3D& v, double r, double h)
{
	// Цилиндр
	double cd = mcGeometry::getDistanceToInfiniteCylinderInside(p, v, r);
	// Плоскости
	double vz = v.z();
	double pd = (vz < 0) ? -p.z() / vz : (vz > 0) ? (h - p.z()) / vz : DBL_MAX;
	return MIN(pd, cd);
}

double mcGeometry::getDistanceToCylinderOutside(const geomVector3D& p, const geomVector3D& v, double r, double h)
{
	double z = p.z(), vz = v.z();
	double R2 = r*r;
	double rr = p.sqLengthXY();

	if (rr <= R2 && (z < 0 || z > h)) // означает, что в направлении Z точка перед или за торцами цилиндра
	{
		// Из-за ошибок округления есть проблема идентификации улетания 
		// на основе нахождения за пределами по Z
		if (z <= h / 2 && vz <= 0) return DBL_MAX;
		else if (z >= h / 2 && vz >= 0) return DBL_MAX;
		else
		{
			double pd = (vz > 0) ? -z / vz : (h - z) / vz;
			geomVector3D pp = p + (v * pd);
			if (pp.sqLengthXY() < R2)
				return pd;
			else
				return DBL_MAX;
		}
	}
	else
	{
		// Пересечение с цилиндром
		double cd = mcGeometry::getDistanceToInfiniteCylinderOutside(p, v, r);
		if (cd == DBL_MAX) return DBL_MAX;

		if (z <= 0) {
			if (fabs(vz) <= DBL_EPSILON) return DBL_MAX;
			double pd = -z / vz;
			geomVector3D pp = p + (v * pd);
			if (pp.sqLengthXY() < R2)
				return pd;
		}
		else if (z >= h)
		{
			if (fabs(vz) <= DBL_EPSILON) return DBL_MAX;
			double pd = (h - z) / vz;
			geomVector3D pp = p + (v * pd);
			if (pp.sqLengthXY() < R2)
				return pd;
		}

		geomVector3D pp = p + (v * cd);
		if (pp.z() <= 0 || pp.z() >= h)
			return DBL_MAX;
		else
			return cd;
	}
}

double mcGeometry::getDistanceToPrismInside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h)
{
	double dx = (v.x() < 0) ? fabs((p.x() + 0.5*ax) / v.x()) : (v.x() > 0) ? fabs(((p.x() - 0.5*ax)) / v.x()) : DBL_MAX;
	double dy = (v.y() < 0) ? fabs((p.y() + 0.5*ay) / v.y()) : (v.y() > 0) ? fabs(((p.y() - 0.5*ay)) / v.y()) : DBL_MAX;
	double dz = (v.z() < 0) ? fabs(p.z() / v.z()) : (v.z() > 0) ? fabs((h - p.z()) / v.z()) : DBL_MAX;
	return MIN(dx, MIN(dy, dz));
}

double mcGeometry::getDistanceToPrismOutside(const geomVector3D& p, const geomVector3D& v, double ax, double ay, double h)
{
	// Упрощаем задачу так, чтобы частица двигалась по всем координатам 
	// в положительном направлении и стартовая координата задавалась
	// относительно угла с наименьшими координатами
	geomVector3D pp(p.x() + 0.5*ax, p.y() + 0.5*ay, p.z());
	geomVector3D vv = v;
	if (vv.x() < 0) { vv(0) = -vv.x(); pp(0) = ax - pp.x(); }
	if (vv.y() < 0) { vv(1) = -vv.y(); pp(1) = ay - pp.y(); }
	if (vv.z() < 0) { vv(2) = -vv.z(); pp(2) = h - pp.z(); }

	// Необходимо проверить только 3 грани на предмет попадания
	if (pp.z() <= 0) {
		if (vv.z() == 0) return DBL_MAX;
		double d = fabs(pp.z() / vv.z());
		geomVector3D c = pp + (vv*d);
		if (c.x() >= 0 && c.x() < ax && c.y() >= 0 && c.y() < ay)
			return d; // попали во фронтальную грань
	}

	if (pp.x() <= 0) {
		if (vv.x() == 0) return DBL_MAX;
		double d = fabs(pp.x() / vv.x());
		geomVector3D c = pp + (vv*d);
		if (c.z() >= 0 && c.z() < h && c.y() >= 0 && c.y() < ay)
			return d;
	}

	if (pp.y() <= 0) {
		if (vv.y() == 0) return DBL_MAX;
		double d = fabs(pp.y() / vv.y());
		geomVector3D c = pp + (vv*d);
		if (c.z() >= 0 && c.z() < h && c.x() >= 0 && c.x() < ax)
			return d;
	}

	// Не попали ни в одну из граней
	return DBL_MAX;

}

double mcGeometry::getDistanceToConeInside(const geomVector3D& p, const geomVector3D& v, double r, double f)
{
	// Решение квадратного уравнения (ad + 2bd + c = 0) 
	// относительно параметра расстояния вдоль линии траектории частицы.

	// Однако, ситуация сложнее простого квадратного уравннения.
	// Он за решение принимает и пересечения с зеркальным отражением конуса.
	// Поэтому, разделяем задачу на отдельные варианты.

	double x = p.x(), y = p.y(), z = p.z();
	double vx = v.x(), vy = v.y(), vz = v.z();

	// Вариант 1. Частица вылетает из конуса в бесконечность.
	double rr = sqrt(r * r + f * f);
	if (vz < 0 && -vz >= f / rr - FLT_EPSILON)
		return DBL_MAX;

	// Параметры квадратного уравнения
	double t = vz * r / f;
	double a = vx * vx + vy * vy - t * t;
	t = 1 - z / f;
	double b = x * vx + y * vy + r * r * t * vz / f;
	double c = x * x + y * y - r * r * t * t;
	double det = b * b - a * c;

	// Вариант 2. a = 0 сооствествует движению вдоль поверхности.
	if (fabs(a) <= FLT_EPSILON)
	{
		if (vz > 0)
			return -0.5 * c / b;
		else
			return DBL_MAX;
	}

	// На всякий случай так как для точек внутри этого не должно быть
	if (det <= 0)
		return DBL_MAX;
	else
		// Если все правильно и частица находится внутри, то единственное правильное
		// большее значение, так как минимальное должно быть отрицатильным.
		// Если минимальное всетаки положительно, то, опять же в нормальном случае,
		// имеет место старт из окрестности поверхности внутрь.
		// Если параметр a оказывается меньше 0, то это означает второе пересечение с зеркальной
		// частью конуса. В этом случае правильным то же будет следующий корень.
		return (-b + sqrt(det)) / a;
}

double mcGeometry::getDistanceToConeOutside(const geomVector3D& p, const geomVector3D& v, double r, double f)
{
	// Сразу отсечь частицы, летящие от конуса
	double s = r / p.lengthXY();
	geomVector3D nn(p.x() * s, p.y() * s, -f);	// нормаль к поверхности конуса в плоскости сечения, проходящей через ось Z и точку частицы
	geomVector3D vv(-nn.y(), nn.x(), 0);
	nn = nn ^ vv;
	s = nn * v;
	if (s >= 0) return DBL_MAX;

	// Параметры квадратного уравнения
	double x = p.x(), y = p.y(), z = p.z();
	double vx = v.x(), vy = v.y(), vz = v.z();
	double t = vz * r / f;
	double a = vx * vx + vy * vy - t * t;
	t = 1 - z / f;
	double b = x * vx + y * vy + r * r * t * vz / f;
	double c = x * x + y * y - r * r * t * t;
	double det = b * b - a * c;

	// Пересечений может не быть вообще
	if (det <= 0) return DBL_MAX;

	double cd;
	if (fabs(a) <= FLT_EPSILON)
		cd = -0.5 * c / b;
	else
	{
		double sd = sqrt(det);
		cd = (-b - sd) / a;
		if (cd < 0)
			cd = (-b + sd) / a;
	}

	// Приемлемо только пересечение с реальным конусом
	geomVector3D cp = p + (v * cd);
	if (cp.z() < f)
		return fabs(cd);
	else
		return DBL_MAX;
}

double mcGeometry::getDistanceToConeSlabInside(const geomVector3D& p, const geomVector3D& v, double z1, double z2, double r1, double r2)
{
	double dd = 0;
	if(r2 == r1) 
		dd = mcGeometry::getDistanceToInfiniteCylinderInside(p, v, r1);
	else
	{
		double x = p.x(), y = p.y(), z = p.z();
		double vx = v.x(), vy = v.y(), vz = v.z();
		double f = (r2 - r1) / (z2 - z1);
		double rfz = r1 - f * z1 + f*z;
		double a = 1 - (1 + f * f) * vz * vz;
		double b = x * vx + y * vy - rfz * f * vz ;
		double c = x * x + y * y - rfz * rfz;

		double det = b * b - a * c;
		if (det <= 0) return DBL_MAX;

		if (fabs(a) <= FLT_EPSILON)
			dd = -0.5 * c / b;
		else
		{
			double sd = sqrt(det);
			double d1, d2;
			if (a > 0) { d1 = (-b - sd) / a; d2 = (-b + sd) / a; }
			else { d2 = (-b - sd) / a; d1 = (-b + sd) / a; }
			dd = (d2 <= 0) ? DBL_MAX : (d1 > 0) ? d1 : d2;
		}
	}
	if(dd == DBL_MAX) return DBL_MAX;
	double z = p.z() + v.z() * dd;
	if (z <= z2 && z >= z1) return dd;
	else return DBL_MAX;
}

double mcGeometry::getDistanceToConeSlabOutside(const geomVector3D& p, const geomVector3D& v, double z1, double z2, double r1, double r2)
{
	double dd = 0;
	if(r2 == r1) 
		dd = mcGeometry::getDistanceToInfiniteCylinderOutside(p, v, r1);
	else
	{
		double x = p.x(), y = p.y(), z = p.z();
		double vx = v.x(), vy = v.y(), vz = v.z();
		double f = (r2 - r1) / (z2 - z1);
		double rfz = r1 - f * z1 + f*z;
		double a = 1 - (1 + f * f) * vz * vz;
		double b = x * vx + y * vy - rfz * f * vz;
		double c = x * x + y * y - rfz * rfz;

		double det = b * b - a * c;
		if (det <= 0) return DBL_MAX;

		if (fabs(a) <= FLT_EPSILON)
			dd = -0.5 * c / b;
		else
		{
			double sd = sqrt(det);
			double d1, d2;
			if (a > 0) { d1 = (-b - sd) / a; d2 = (-b + sd) / a; }
			else { d2 = (-b - sd) / a; d1 = (-b + sd) / a; }
			dd = (d2 <= 0) ? DBL_MAX : (d1 > 0) ? d1 : d2;
		}
	}
	if (dd == DBL_MAX) return DBL_MAX;
	double z = p.z() + v.z() * dd;
	if (z <= z2 && z >= z1) return dd;
	else return DBL_MAX;
}

double mcGeometry::getDistanceToRectanglePipeInside(const geomVector3D& p, const geomVector3D& v,
	double x1, double x2, double y1, double y2)
{
	double dx1 = p.x() - x1;
	double cx1 = v.x() >= 0 ? DBL_MAX : dx1 <= 0 ? 0 : -dx1 / v.x();

	double dx2 = x2 - p.x();
	double cx2 = v.x() <= 0 ? DBL_MAX : dx2 <= 0 ? 0 : dx2 / v.x();

	double dy1 = p.y() - y1;
	double cy1 = v.y() >= 0 ? DBL_MAX : dy1 <= 0 ? 0 : -dy1 / v.y();

	double dy2 = y2 - p.y();
	double cy2 = v.y() <= 0 ? DBL_MAX : dy2 <= 0 ? 0 : dy2 / v.y();

	return MIN(MIN(cx1, cx2), MIN(cy1, cy2));
}

double mcGeometry::getDistanceToRectanglePipeOutside(const geomVector3D& p, const geomVector3D& v,
	double x1, double x2, double y1, double y2)
{
	// Упрощаем задачу так, чтобы частица двигалась по всем координатам 
	// в положительном направлении и стартовая координата задавалась
	// относительно угла с наименьшими координатами
	geomVector3D pp(p.x() - x1, p.y() - y1, p.z());
	geomVector3D vv = v;
	if (vv.x() < 0) { vv(0) = -vv.x(); pp(0) = x2 - x1 - pp.x(); }
	if (vv.y() < 0) { vv(1) = -vv.y(); pp(1) = y2 - y1 - pp.y(); }

	// Необходимо проверить только 2 грани на предмет попадания
	if (pp.x() <= 0) {
		if (vv.x() == 0) return DBL_MAX;
		double d = fabs(pp.x() / vv.x());
		geomVector3D c = pp + (vv*d);
		if (c.y() >= 0 && c.y() < y2 - y1)
			return d;
	}

	if (pp.y() <= 0) {
		if (vv.y() == 0) return DBL_MAX;
		double d = fabs(pp.y() / vv.y());
		geomVector3D c = pp + (vv*d);
		if (c.x() >= 0 && c.x() < x2 - x1)
			return d;
	}

	// Не попали ни в одну из граней
	return DBL_MAX;
}

double mcGeometry::getDistanceToRectangleConeInside(const geomVector3D& p, const geomVector3D& v,
	double x, double cosx, double sinx,
	double y, double cosy, double siny, double z)
{
	double px = p.x(), py = p.y(), pz = p.z();
	double vx = v.x(), vy = v.y(), vz = v.z();

	// Используем симметричность и позиционируем частицу в левый нижний квадрант плоскости XY
	if (vx < 0) { px = -px; vx = -vx; }
	if (vy < 0) { py = -py; vy = -vy; }

	// Скалярное произведение вектора скорости и нормали к грани
	double vn = vx * cosx + vz * sinx;
	double cx1 = vn >= 0 ? DBL_MAX : (-(px + x) * cosx - (pz - z) * sinx) / vn;
	vn = -vx * cosx + vz * sinx;
	double cx2 = vn >= 0 ? DBL_MAX : ((px - x) * cosx - (pz - z) * sinx) / vn;

	vn = vy * cosy + vz * siny;
	double cy1 = vn >= 0 ? DBL_MAX : (-(py + y) * cosy - (pz - z) * siny) / vn;
	vn = -vy * cosy + vz * siny;
	double cy2 = vn >= 0 ? DBL_MAX : ((py - y) * cosy - (pz - z) * siny) / vn;

	return MIN(MIN(cx1, cy1), MIN(cx2, cy2));
}

double mcGeometry::getDistanceToRectangleConeOutside(const geomVector3D& p, const geomVector3D& v,
	double x, double cosx, double sinx,
	double y, double cosy, double siny, double z)
{
	double px = p.x(), py = p.y(), pz = p.z();
	double vx = v.x(), vy = v.y(), vz = v.z();

	// Используем симметричность и позиционируем частицу в левый нижний квадрант плоскости XY
	if (px > 0) { px = -px; vx = -vx; }
	if (py > 0) { py = -py; vy = -vy; }

	double vn = -vx * cosx - vz * sinx;
	double cx = vn >= 0 ? -1 : ((px + x) * cosx + (pz - z) * sinx) / vn;

	vn = -vy * cosy - vz * siny;
	double cy = vn >= 0 ? -1 : ((py + y) * cosy + (pz - z) * siny) / vn;

	if (cx >= 0)
	{
		double dz = pz + vz * cx - z;
		if (fabs(py + vy * cx) <= (y + dz * siny / cosy))
			return cx;
	}

	if (cy >= 0)
	{
		double dz = pz + vz * cy - z;
		if (fabs(px + vx * cy) <= (x + dz * sinx / cosx))
			return cy;
	}

	return DBL_MAX;
}

double mcGeometry::TrLenInVoxel(double x1, double y1, double z1
	, double x2, double y2, double z2
	, double xv1, double yv1, double zv1
	, double xv2, double yv2, double zv2)
{
	double f;
	// Make partical direction positive along Z-axis.
	if (z2 < z1) {
		f = x1; x1 = x2; x2 = f;
		f = y1; y1 = y2; y2 = f;
		f = z1; z1 = z2; z2 = f;
	}
	//Find partical travel limits between planes zv1 < z < zv2 
	if (z1 >= zv2 || z2 <= zv1)
		return 0;
	if (z1 < zv1) {
		f = (zv1 - z1) / (z2 - z1);
		x1 = x1 + (x2 - x1)*f;
		y1 = y1 + (y2 - y1)*f;
		z1 = zv1;
	}
	if (z2 > zv2) {
		f = (zv2 - z1) / (z2 - z1);
		x2 = x1 + (x2 - x1)*f;
		y2 = y1 + (y2 - y1)*f;
		z2 = zv2;
	}
	// Swap X coordinates if necessary
	// The rectangle geometry of voxel allows that
	if (x2 < x1) {
		f = x1; x1 = x2; x2 = f;
	}
	// Find intersections partical travel limits between planes xv1 < x < xv2 
	if (x1 >= xv2 || x2 <= xv1)
		return 0;
	if (x1 < xv1) {
		f = (xv1 - x1) / (x2 - x1);
		z1 = z1 + (z2 - z1)*f;
		y1 = y1 + (y2 - y1)*f;
		x1 = xv1;
	}
	if (x2 > xv2) {
		f = (xv2 - x1) / (x2 - x1);
		z2 = z1 + (z2 - z1)*f;
		y2 = y1 + (y2 - y1)*f;
		x2 = xv2;
	}
	// Swap Y coordinates if necessary
	if (y2 < y1) {
		f = y1; y1 = y2; y2 = f;
	}
	// Find intersections partical travel limits between planes yv1 < y < yv2
	if (y1 >= yv2 || y2 <= yv1)
		return 0;
	if (y1 < yv1) {
		f = (yv1 - y1) / (y2 - y1);
		x1 = x1 + (x2 - x1)*f;
		z1 = z1 + (z2 - z1)*f;
		y1 = yv1;
	}
	if (y2 > yv2) {
		f = (yv2 - y1) / (y2 - y1);
		x2 = x1 + (x2 - x1)*f;
		z2 = z1 + (z2 - z1)*f;
		y2 = yv2;
	}
	// Calculate the final length
	return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
}

double mcGeometry::getDistanceToSphereInside(const geomVector3D& p, const geomVector3D& v, double r)
{
	double b = p * v;
	double c = p * p - r * r;
	double det = b * b - c;
	if (det < 0) det = 0;      // eliminate surface vicinity problem
	// when points is inside, only one positive solution exists
	double dist = -b + sqrt(det);
	return dist >= 0 ? dist : -dist;
}

double mcGeometry::getDistanceToSphereOutside(const geomVector3D& p, const geomVector3D& v, double r)
{
	double b = p * v;
	// if moves out of sphere
	if (b >= 0)
		return DBL_MAX;

	double c = p * p - r * r;
	double det = b*b - c;
	if (det <= 0) // no intersection
		return DBL_MAX;

	// Only two solutions can exist. Smallest distance is what we are looking for.
	// If it is negative, then we move out of sphere.
	double dist = -b - sqrt(det);
	return dist >= 0 ? dist : -dist;
}

double mcGeometry::getDistanceToConvexPolygonCircleInside(const geomVector3D& p, const geomVector3D& v, const std::vector<double>& pz, const std::vector<double>& pr)
{
	// Идем по слоям Z.
	// Если пересечение в пределах радиуса очередного слоя, то накапливаем и двигаемся дальше.
	// Если нет, то очередное пересечение будет с конусом.

	unsigned i;
	double dcount = 0;
	geomVector3D pcurrent(p);
	double z = p.z();
	for (i = 1; i < pz.size(); i++) if (z <= pz[i] + DBL_EPSILON) break;
	double z1 = pz[i - 1], z2 = pz[i];
	double r1 = pr[i - 1], r2 = pr[i];

	double vz = v.z();
	// Движение параллельно плоскости XY
	if (vz == 0)
	{
		double r = r1 + (r2 - r1) * (z - z1) / (z2 - z1);
		double b = (p.x()*v.x() + p.y()*v.y());
		double c = (p.sqLengthXY() - r*r);
		double det = b*b - c;
		double cd = 0;
		if (det > 0)
			cd = -b + sqrt(det);  // when points is inside, only one positive solution exists
		return (cd < 0 ? 0 : cd);
	}

	while (true)
	{
		if (vz < 0)
		{
			double dd = (z1-z) / vz;
			geomVector3D pc = pcurrent + (v * dd);
			double rr = pc.lengthXY();
			if (rr <= r1)
			{
				dcount += dd;
				if (i <= 1)
					break;
				i--;
				pcurrent = pc;
				z = z1;
				z2 = z1;
				r2 = r1;
				z1 = pz[i - 1];
				r1 = pr[i - 1];
				continue;
			}
		}
		else
		{
			double dd = (z2 - z) / vz;
			geomVector3D pc = pcurrent + (v * dd);
			double rr = pc.lengthXY();
			if (rr <= r2)
			{
				dcount += dd;
				if (i >= pr.size() - 1)
					break;
				i++;
				pcurrent = pc;
				z = z2;
				z1 = z2;
				r1 = r2;
				z2 = pz[i];
				r2 = pr[i];
				continue;
			}
		}

		// К этому моменту остается только конус
		double dd = getDistanceToConeSlabInside(pcurrent, v, z1, z2, r1, r2);
		// Hack!!! Нельзя возвращать бесконечность внутри.
		// Вероятно причина в поверхностных эффектах.
		if (dd != DBL_MAX) 
			dcount += dd;
		break;
	}

	return dcount;
}

double mcGeometry::getDistanceToConvexPolygonCircleOutside(const geomVector3D & p, const geomVector3D & v, const std::vector<double>& pz, const std::vector<double>& pr)
{
	// Аналогично идем по слоям Z.
	// В зависимости от места пересечения определяем пересечение с боковой поверхностью
	// или движемся дальше.

	int i = -1;
	double dcount = 0;
	geomVector3D pcurrent(p);
	double z = p.z();
	double vz = v.z();
	double z1 = pz[0], z2 = pz.back();
	double r1 = pr[0], r2 = pr.back();

	if (z <= z1 + DBL_EPSILON)
	{
		if (vz <= 0)
			return DBL_MAX;

		double dd = (z1 - z) / vz;
		geomVector3D pc = pcurrent + (v * dd);
		if (pc.lengthXY() <= r1)
			return dd;
		else
		{
			dcount += dd;
			i = 1;
			pcurrent = pc;
			z2 = pz[i];
			r2 = pr[i];
		}
	}
	else if (z >= z2 - DBL_EPSILON)
	{
		if (vz >= 0)
			return DBL_MAX;

		double dd = (z2 - z) / vz;
		geomVector3D pc = pcurrent + (v * dd);
		if (pc.lengthXY() <= r2)
			return dd;
		else
		{
			dcount += dd;
			i = int(pz.size() - 1);
			pcurrent = pc;
			z1 = pz[i-1];
			r1 = pr[i-1];
		}
	}

	// Теперь находимся точно между крайними слоями.
	// При необходимости определяемся между какими.
	if (i == -1)
	{
		for (i = 1; i < pz.size(); i++) if (z <= pz[i] + DBL_EPSILON) break;
		z1 = pz[i - 1]; z2 = pz[i];
		r1 = pr[i - 1]; r2 = pr[i];
	}

	// Движение параллельно плоскости XY
	if (vz == 0)
	{
		double r = r1 + (r2 - r1) * (z - z1) / (z2 - z1);
		double b = (p.x()*v.x() + p.y()*v.y());
		double c = (p.sqLengthXY() - r*r);
		double det = b*b - c;
		if (det <= 0) return DBL_MAX;
		double cd = (-b - sqrt(det));
		return (cd < 0) ? DBL_MAX : cd;
	}

	while (true)
	{
		double dd = getDistanceToConeSlabOutside(pcurrent, v, z1, z2, r1, r2);
		if (dd == DBL_MAX)
		{
			if (vz < 0)
			{
				i--;
				if (i <= 0) return DBL_MAX;
				dd = (z1 - z) / vz;
				pcurrent += (v * dd);
				z = z1;
				z2 = z1; r2 = r1;
				z1 = pz[i - 1];
				r1 = pr[i - 1];
				dcount += dd;
			}
			else
			{
				i++;
				if (i >= pz.size()) return DBL_MAX;
				dd = (z2 - z) / vz;
				pcurrent += (v * dd);
				z = z2;
				z1 = z2; r1 = r2;
				z2 = pz[i];
				r2 = pr[i];
				dcount += dd;
			}
		}
		else {
			dcount += dd; break;
		}
	}

	return dcount;
}
