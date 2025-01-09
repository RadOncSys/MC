#include "mcTransportMantleBlock.h"
#include "mcGeometry.h"
#include "mcDefs.h"
#include "../Geometry/vec2d.h"
#include <float.h>

mcTransportMantleBlock::mcTransportMantleBlock(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, 
	double r1, double h, std::vector<double>& plgnX, std::vector<double>& plgnY)
	: mcTransport(orgn, z, x), x_(plgnX), y_(plgnY), r1_(r1), h_(h)
{
	nsides_ = (int)plgnX.size();
	geomVector3D vy(0, 0, 1);
	sides_ = std::make_unique<std::vector<mcGeomRectSide>>(nsides_);

	// У полигона должна быть ориентация против часовой стрелки.
	// Берем крайнюю левую точку и смотрим куда идет следующая.
	double xmin = DBL_MAX;
	int imin = 0;
	for (int i = 0; i < nsides_; i++)
	{
		if (x_[i] < xmin)
		{
			xmin = x_[i];
			imin = i;
		}
	}
	if (y_[(imin + nsides_ - 1) % nsides_] < y_[imin] ||
		y_[(imin + 1) % nsides_] > y_[imin])
	{
		// Переориентация
		for (int i = 0; i < nsides_ / 2; i++)
		{
			double t = x_[i];
			x_[i] = x_[nsides_ - i - 1];
			x_[nsides_ - i - 1] = t;
			t = y_[i];
			y_[i] = y_[nsides_ - i - 1];
			y_[nsides_ - i - 1] = t;
		}
	}

	for (int i = 0; i < nsides_; i++)
	{
		int i1 = (i + 1) % nsides_;
		double dx = x_[i1] - x_[i];
		double dy = y_[i1] - y_[i];
		double ax = sqrt(dx * dx + dy * dy);
		geomVector3D vx(dx, dy, 0);
		geomVector3D p((x_[i1] + x_[i]) / 2, (y_[i1] + y_[i]) / 2, h_ / 2);
		sides_->at(i).SetGeometry(p, vx, vy, ax, h_);
	}
}

double mcTransportMantleBlock::getDistanceInside(mcParticle& p)	const
{
	// Внешний цилиндр
	double cd1 = mcGeometry::getDistanceToInfiniteCylinderInside(p.p, p.u, r1_);

	// Внутренний полигональный цилиндр
	double cd2 = DBL_MAX;
	for (int i = 0; i < nsides_; i++)
	{
		double d = sides_->at(i).getDistance(p.p, p.u);
		if (d < cd2) cd2 = d;
	}

	// Торцевые плоскости
	double vz = p.u.z();
	double pd = (vz < 0) ? -p.p.z() / vz : (vz > 0) ? (h_ - p.p.z()) / vz : DBL_MAX;

	return MIN(MIN(cd1, cd2), pd);
}

double mcTransportMantleBlock::getDistanceOutside(mcParticle& p) const
{
	// Отсекаем движение от торцевых плоскостей
	double pz = p.p.z(), uz = p.u.z();
	if ((pz < 0 && uz <= 0) || (pz > h_ && uz >= 0))
		return DBL_MAX;

	// Отсекаем столкновения с боковой поверхностью цилинда
	if (p.p.lengthXY() > r1_)
	{
		double cd = mcGeometry::getDistanceToInfiniteCylinderOutside(p.p, p.u, r1_);
		if (cd == DBL_MAX)
			return DBL_MAX;
		else
		{
			auto pp = p.p + (p.u * cd);
			if (pp.z() <= h_ && pp.z() >= 0)
				return cd;
			else if (pz >= 0 && pz <= h_)
				return DBL_MAX;
		}
	}

	// Текущий счетчик перемещений и положения частицы
	double dist = 0;
	auto c = p.p;

	// Проверяем попдание в торцы цилиндра и в случае успеха 
	// перемещаем частицу на его поверхность
	if (c.z() <= 0)
	{
		dist = -pz / uz;
		c += p.u * (dist + MINDELTA);
		if (c.lengthXY() > r1_)
			return DBL_MAX;
	}
	else if (c.z() >= h_)
	{
		dist = -(pz - h_) / uz;
		c += p.u * (dist + MINDELTA);
		if (c.lengthXY() > r1_)
			return DBL_MAX;
	}

	// Если мы переместились на поверхность, то отсекаем ситуацию что мы уже в теле объекта
	if (dist > 0 && !isPointInPlgn(c.x(), c.y()))
		return dist;
	else
	{
		double cd2 = DBL_MAX;
		if (c.z() >= 0 && c.z() <= h_)
		{
			for (int i = 0; i < nsides_; i++)
			{
				double d = sides_->at(i).getDistance(c, p.u);
				if (d < cd2) cd2 = d;
			}
		}
		return cd2 == DBL_MAX ? DBL_MAX : cd2 + dist;
	}

	// Этого выхода быть не должно.
	// Если добрались, то есть какая-то проблема.
	cout << "Something wrong in mcTransportMantleBlock::getDistanceOutside" << endl;

	return DBL_MAX;
}

bool mcTransportMantleBlock::isPointInPlgn(double x, double y) const
{
	// Собираем все сечения полигона линиями y = const, где y - координата тестируемой точки.

	// Счетчик пересечений слева от точки
	int count = 0;

	// Параметры состояния предыдущего совпадения точки полигона с линией y
	double xprev = 0;
	int dir = 1;

	// Стартовать нужно с с хорошего сегмента
	int i0 = 0;
	for (int idx = 0; idx < nsides_ - 1; idx++)
	{
		if ((y - y_[idx]) * (y - y_[idx + 1]) < 0)
		{
			i0 = idx;
			break;
		}
	}

	for (int idx = 0; idx < nsides_; idx++)
	{
		int i = (i0 + idx) % nsides_;
		int i1 = (i0 + idx + 1) % nsides_;
		double dy = y - y_[i];
		double dy1 = y - y_[i1];
		if (dy * dy1 > 0) continue;
		else if (dy * dy1 < 0)
		{
			double cx = x_[i] + (x_[i1] - x_[i]) * (y - y_[i]) / (y_[i1] - y_[i]);
			if (cx < x)
				count++;
		}
		else if (dy == 0 && dy1 == 0) continue;
		else if (dy == 0)
		{
			if (dy1 * dir > 0)
			{
				if (x >= (xprev + x_[i]) / 2) 
					count++;
			}
		}
		else if (dy1 == 0)
		{
			xprev = x_[i];
			dir = dy > 0 ? 1 : -1;
		}
	}

	return (count % 2) > 0;
}

double mcTransportMantleBlock::getDNearInside(const geomVector3D& p) const
{
	return 0;
}

void mcTransportMantleBlock::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Radius:\t" << r1_ << endl;
	os << "Height:\t" << h_ << endl;
	os << "X: ";
	for (int i = 0; i < nsides_; i++)
		os << "\t" << x_[i];
	os << endl;
	os << "Y: ";
	for (int i = 0; i < nsides_; i++)
		os << "\t" << y_[i];
	os << endl;
}

void mcTransportMantleBlock::dumpVRML(ostream& os) const
{
	// Единственная проблема - это заполнение торцов.
	// Идея решения в следующем.
	// Боковые сторны заполняем как обычно, но не используем готовые функции рисования цилиндров.
	// Это нужно для контроля над точками полигона круга, которые должны быть согласованы с торцами.
	// Находим центр внутреннего полигона из которого будем проводить линии для определения 
	// какие точки двух полигонов нужно соединять для построения треугольников торцов.

	// Полигон внешнего круга.
	int i, na = 72;
	double da = 2 * PI / na;
	std::vector<double> rx(na);
	std::vector<double> ry(na);
	for (i = 0; i < na; i++)
	{
		rx[i] = r1_ * cos(i * da);
		ry[i] = r1_ * sin(i * da);
	}

	// Центр внутреннего полигона
	double x1 = DBL_MAX, x2 = -DBL_MAX;
	double y1 = DBL_MAX, y2 = -DBL_MAX;
	for (i = 0; i < nsides_; i++)
	{
		if (x_[i] < x1) x1 = x_[i];
		if (x_[i] > x2) x2 = x_[i];
		if (y_[i] < y1) y1 = y_[i];
		if (y_[i] > y2) y2 = y_[i];
	}
	geomVector2D P0((x1 + x2) / 2,  (y1 + y2) / 2);

	os << "# Mantle Block: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	// Боковая стенка цилиндра

	os << "    Transform {" << endl;
	os << "      children Shape {" << endl;
	os << "        appearance Appearance {" << endl;
	os << "          material Material {" << endl;
	os << "            diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "            transparency " << transparancy_ << endl;
	os << "          }" << endl;
	os << "        }" << endl;
	os << "        geometry IndexedFaceSet {" << endl;
	os << "            coord Coordinate {" << endl;
	os << "                point [" << endl;

	for (i = 0; i < na; i++) {
		geomVector3D p = geomVector3D(rx[i], ry[i], 0) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;
		p = geomVector3D(rx[i], ry[i], h_) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		if (i < na - 1) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	for (i = 0; i < na; i++) {
		os << "                " << 2 * i << ", " << 2 * ((i + 1) % na) << ", " << 2 * ((i + 1) % na) + 1 << ", " << 2 * i + 1;
		if (i < na - 1) os << ", -1,";
		os << endl;
	}

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;

	// Внутренняя стенка

	os << "    Transform {" << endl;
	os << "      children Shape {" << endl;
	os << "        appearance Appearance {" << endl;
	os << "          material Material {" << endl;
	os << "            diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "            transparency " << transparancy_ << endl;
	os << "          }" << endl;
	os << "        }" << endl;
	os << "        geometry IndexedFaceSet {" << endl;
	os << "            coord Coordinate {" << endl;
	os << "                point [" << endl;

	for (i = 0; i < nsides_; i++) {
		geomVector3D p = geomVector3D(x_[i], y_[i], 0) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;
		p = geomVector3D(x_[i], y_[i], h_) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		if (i < nsides_ - 1) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	for (i = 0; i < na; i++) {
		os << "                " << 2 * i << ", " << 2 * i + 1 << ", " << 2 * ((i + 1) % na) + 1 << ", " << 2 * ((i + 1) % na);
		if (i < nsides_ - 1) os << ", -1,";
		os << endl;
	}

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;

	// Нижний и верхний торец сразу. Точки очевидны. Магия в индексах

	os << "    Transform {" << endl;
	os << "      children Shape {" << endl;
	os << "        appearance Appearance {" << endl;
	os << "          material Material {" << endl;
	os << "            diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "            transparency " << transparancy_ << endl;
	os << "          }" << endl;
	os << "        }" << endl;
	os << "        geometry IndexedFaceSet {" << endl;
	os << "            coord Coordinate {" << endl;
	os << "                point [" << endl;

	for (i = 0; i < na; i++) {
		geomVector3D p = geomVector3D(rx[i], ry[i], 0) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;
		p = geomVector3D(rx[i], ry[i], h_) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		//if (i < na - 1) os << ", ";
		os << endl;
	}

	for (i = 0; i < nsides_; i++) {
		geomVector3D p = geomVector3D(x_[i], y_[i], 0) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;
		p = geomVector3D(x_[i], y_[i], h_) * mttow_;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		if (i < na - 1) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	// Собственно магия. 
	// Проводим линии от центра внутреннего полигона последовательно в точки круга и 
	// отслеживаем индексы и переходы узлов внутреннего полигона.
	// В зависимости от состояния определяем треугольники.

	int ir_current = na - 1;

	// Первая точка полигона слева от направления на первую точку круга
	int ip_current = 0;
	geomVector2D pr0(rx[ir_current], ry[ir_current]);
	geomVector2D vr0 = pr0 - P0;
	geomVector2D nr0 = vr0;
	nr0.turnLeft();

	for (i = 0; i < nsides_; i++)
	{
		int i1 = (i + 1) % nsides_;
		geomVector2D vp = geomVector2D(x_[i], y_[i]) - P0;
		geomVector2D vp1 = geomVector2D(x_[i1], y_[i1]) - P0;
		if ((vp * nr0) * (vp1 * nr0) <= 0 && (vp * pr0) > 0)
		{
			ip_current = i1;
			break;
		}
	}

	for (int i = 0; i < na; i++) 
	{
		// Направление на очередную точку круга
		pr0.set(rx[ir_current], ry[ir_current]);
		vr0 = pr0 - P0;
		nr0 = vr0;
		nr0.turnLeft();

		// Проверяем, не пересекли ли границу сегмента внутреннего полигона
		bool isCrossed = false;
		for (int j = 0; j < nsides_; i++)
		{
			geomVector2D vp = geomVector2D(x_[ip_current], y_[ip_current]) - P0;
			if((vp * nr0) > 0)
				break;
			isCrossed = true;
			int ip_next = (ip_current + 1) % nsides_;

			os << "                " << 2 * ir_current << ", " << 2 * (na + ip_current) << ", " << 2 * (na + ip_next);
			os << ", -1," << endl;
			os << "                " << 2 * ir_current + 1 << ", " << 2 * (na + ip_next) + 1 << ", " << 2 * (na + ip_current) + 1;
			os << ", -1," << endl;

			ip_current = ip_next;
		}

		// Независимо от того пересекли ли сегмент полигона или нет добавляем треугольник с основанием на круге.
		os << "                " << 2 * i << ", " << 2 * ir_current << ", " << 2 * (na + ip_current);
		os << ", -1," << endl;
		os << "                " << 2 * ir_current + 1 << ", " << 2 * i + 1 << ", " << 2 * (na + ip_current) + 1;
		if (i < na - 1) os << ", -1,";
		os << endl;

		ir_current = i;
	}

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;

	os << "  ]" << endl;
	os << "}" << endl;
}
