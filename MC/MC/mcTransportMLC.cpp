#include "mcTransportMLC.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportMLC::mcTransportMLC(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
	double focus, double r, double h, double w, double l, double fsx1, double fsx2, double fsy1, double fsy2)
	: mcTransport(orgn, vz, vx)
	, focus_(focus), r_(r), h_(h), w_(w), l_(l), fsx1_(fsx1), fsx2_(fsx2), fsy1_(fsy1), fsy2_(fsy2)
	, pymin_(0, -w / 2, 0), pymax_(0, w / 2, 0), py1_(0, fsy1, 0), py2_(0, fsy2, 0)
{
	r2_ = r_ * r_;
	dx_ = r_ * (1 - sqrt(1.0 - 0.25 * h_ * h_ / r2_));
	sfy_ = focus_ / (focus_ - h_);

	zxn1_.set(0, focus_, -w_ / 2);
	zxn2_.set(0, -focus_, -w_ / 2);
	zxn1_.normalize();
	zxn2_.normalize();
}

mcTransportMLC::~mcTransportMLC(void)
{
}

double mcTransportMLC::getDistanceInside(mcParticle& p) const
{
	// Идея решения в разбиении на три слоя (шесть кусков) и расчеты по отдельности в каждом

	geomVector3D position(p.p);
	geomVector3D direction(p.u);
	double x = position.x();
	double y = position.y();
	double z = position.z();
	double xx = fsx1_;
	bool isRight = false;
	double dist = 0;

	// Для определения текущего блока проектируем положение частицы на нижнюю плоскость
	double y0 = y * focus_ / (focus_ - z);

	// Определяем содержащий блок и необходимость перевернуть стороны.
	if (y0 < fsy1_ || y0 > fsy2_)
	{
		if (x > 0) isRight = true;
	}
	else
	{
		if (x > 0.5 * (fsx1_ + fsx2_)) isRight = true;
	}

	// Перевернуть если необходимо
	if (isRight)
	{
		position(0) = -position.x();
		direction(0) = -direction.x();
		xx = -fsx2_;
	}

	// Определение растояний в зависимости от блока
	if (y0 < fsy1_)
	{
		dist = getDistanceToBlockInside(position, direction, -0.5*w_, fsy1_, 0);
		position += direction * dist;

		// Если оказались на поверхности следующего блока, то двигаемся еще и в нем
		if (isOnOnTheLeftLeaf(position, fsy1_, xx))
		{
			double d = getDistanceToBlockInside(position, direction, fsy1_, fsy2_, xx);
			position += direction * d;
			dist += d;

			if (isOnOnTheLeftLeaf(position, fsy2_, 0))
				dist += getDistanceToBlockInside(position, direction, fsy2_, 0.5*w_, 0);
		}
	}
	else if (y0 < fsy2_)
	{
		dist = getDistanceToBlockInside(position, direction, fsy1_, fsy2_, xx);
		position += direction * dist;

		if (isOnOnTheLeftLeaf(position, fsy2_, 0))
			dist += getDistanceToBlockInside(position, direction, fsy2_, 0.5*w_, 0);

		else if (isOnOnTheLeftLeaf(position, fsy1_, 0))
			dist += getDistanceToBlockInside(position, direction, -0.5*w_, fsy1_, 0);
	}
	else
	{
		dist = getDistanceToBlockInside(position, direction, fsy2_, 0.5*w_, 0);
		position += direction * dist;

		if (isOnOnTheLeftLeaf(position, fsy2_, xx))
		{
			double d = getDistanceToBlockInside(position, direction, fsy1_, fsy2_, xx);
			position += direction * d;
			dist += d;

			if (isOnOnTheLeftLeaf(position, fsy1_, 0))
				dist += getDistanceToBlockInside(position, direction, -0.5*w_, fsy1_, 0);
		}
	}

	return dist;
}

double mcTransportMLC::getDistanceOutside(mcParticle& p) const
{
	geomVector3D position(p.p);
	geomVector3D direction(p.u);
	double z = position.z();
	double vz = direction.z();
	double dist = 0;

	// Если мы не между слоями Z, то перемещаемся на первый на пути
	if (z >= h_)
	{
		if (vz < 0)
		{
			dist = (h_ - z) / vz;
			position += (direction * dist);
			if (isInZPlaneHit(position.x(), position.y() * sfy_))
				return dist;
		}
		else
			return DBL_MAX;
	}
	else if (z <= 0)
	{
		if (vz > 0)
		{
			dist = -z / vz;
			position += (direction * dist);
			if (isInZPlaneHit(position.x(), position.y()))
				return dist;
		}
		else
			return DBL_MAX;
	}

	// К данному моменту частица находится по Z между 0 и h_
	double dd = isOutsideMLC(position.x(), position.y() * sfy_) ?
		getDistanceToMlcOutsideField(position, direction) :
		getDistanceToMlcInsideField(position, direction);

	return dd == DBL_MAX ? DBL_MAX : dd + dist;
}

double mcTransportMLC::getDNearInside(const geomVector3D& p) const
{
	return 0;
}

double mcTransportMLC::getDistanceToBlockInside(const geomVector3D& p, const geomVector3D& u, double y1, double y2, double xx) const
{
	double cd1 = DBL_MAX, cd2 = DBL_MAX, cd3 = DBL_MAX, cd4 = DBL_MAX, cd5 = DBL_MAX;

	// Расстояние до плоскостей слоя MLC
	double z = p.z(), vz = u.z();
	if (vz < 0) cd1 = -z / vz;
	else if (vz > 0) cd1 = (h_ - z) / vz;

	// Пересечение с наклонными плоскостями XZ
	geomVector3D zxn1(0, focus_, y1), zxn2(0, -focus_, -y2);
	zxn1.normalize();
	zxn2.normalize();
	geomVector3D py1(0, y1, 0), py2(0, y2, 0);

	double dir = -(u * zxn1);
	if (dir > 0)
		cd2 = ((p - py1) * zxn1) / dir;

	dir = -(u * zxn2);
	if (dir > 0)
		cd3 = ((p - py2) * zxn2) / dir;

	// Расстояние до левого торца x = const
	double x = p.x(), vx = u.x();
	if (vx < 0) cd4 = (-l_ + xx - x) / vx;

	// Расстояние до скругленного торца
	// Можно грубо отсечь пересечения за пределами слоя МЛК на основании пересечения с плоскостями z = const
	if (cd1 == DBL_MAX || ((x + vx * cd1) > (xx - dx_)))
	{
		double d = 0;
		double rx = x + r_ - xx;
		z -= h_ / 2;
		double rr = rx * rx + z * z;
		geomVector3D pp(rx, z, 0);
		geomVector3D vv(vx, vz, u.y());
		
		// Потенциальная проблема в том, что можем находиться за пределеами цилиндра.
		// На этот случай проверяем где находимся и при необходимости переезжаем в плоскость через ось цилиндра.
		// Перемещаем на поверхность цилиндра снаружи если необходимо
		if (rr > r2_)
		{
			d = mcGeometry::getDistanceToInfiniteCylinderOutside(pp, vv, r_);
			if (d != DBL_MAX)
				pp += u * d;
		}
		if (d != DBL_MAX)
		{
			double dd = mcGeometry::getDistanceToInfiniteCylinderInside(pp, vv, r_);
			if (dd != DBL_MAX && (rx + dd * vx) >= dx_) cd5 = dd + d;
		}
	}

	double dist = MIN(MIN(cd1, MIN(cd2, cd3)), MIN(cd4, cd5));
	return dist == DBL_MAX ? DBL_MAX : dist;
}

double mcTransportMLC::getDistanceToMlcOutsideField(const geomVector3D& p, const geomVector3D& u) const
{
	// Для начала впихиваем между плоскостями XZ
	double dist = 0;
	double y = p.y() * focus_ / (focus_ - p.z());
	geomVector3D pp(p);

	if (y < -w_ / 2)
	{
		double dir = u * zxn1_;
		if (dir <= 0) return DBL_MAX;
		dist = ((pymin_ - p) * zxn1_) / dir;
		pp += u * dist;
		if (pp.z() >= h_ || pp.z() <= 0)
			return DBL_MAX;
		else if (isOnOutsideLeaf(pp))
			return dist;
		else if (fabs(pp.x()) < dx_)
		{
			// Попали в щель между закрытыми лепестками
			double d = getDistanceToMlcInsideField(p, u);
			return d == DBL_MAX ? DBL_MAX : d;
		}
	}
	else if (y > w_ / 2)
	{
		double dir = u * zxn2_;
		if (dir <= 0) return DBL_MAX;
		dist = ((pymax_ - p) * zxn2_) / dir;
		pp += u * dist;
		if (pp.z() >= h_ || pp.z() <= 0)
			return DBL_MAX;
		else if (isOnOutsideLeaf(pp))
			return dist;
		else if (fabs(pp.x()) < dx_)
		{
			// Попали в щель между закрытыми лепестками
			double d = getDistanceToMlcInsideField(p, u);
			return d == DBL_MAX ? DBL_MAX : d;
		}
	}

	// Теперь остается только одно направление возможных пересечений - ось X
	double x = pp.x();
	if (x < 0)
	{
		// Сначала переводим за плоскость выступающей части
		if (fsx1_ < 0)	// выступает центральная часть 
		{
			if (x <= -l_ + fsx1_)
			{
				if (u.x() <= 0) return DBL_MAX;
				double d = (fsx1_ - l_ - x) / u.x();
				pp += u * d;
				dist += d;
				if (pp.z() <= 0 || pp.z() >= h_) return DBL_MAX;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if (y <= -w_ / 2 && y >= w_ / 2) return DBL_MAX;
				else if (y >= fsy1_ && y <= fsy2_) return dist;
			}
			// Транспорт в двух угловых закутках
			if (y < 0)
			{
				double d = getDistanceToMlcRectangleConeInside(pp, u, -l_, -l_ + fsx1_, -w_ / 2, fsy1_);
				if (d == DBL_MAX) return DBL_MAX;
				pp += u * d;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if (y <= -w_ / 2) return DBL_MAX;
				else return d + dist;
			}
			else
			{
				double d = getDistanceToMlcRectangleConeInside(pp, u, -l_, -l_ + fsx1_, fsy2_, w_ / 2);
				if (d == DBL_MAX) return DBL_MAX;
				pp += u * d;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if (y >= w_ / 2) return DBL_MAX;
				else return d + dist;
			}
		}
		else
		{
			if (x <= -l_)
			{
				if (u.x() <= 0) return DBL_MAX;
				double d = (-l_ - x) / u.x();
				pp += u * d;
				dist += d;
				if (pp.z() <= 0 || pp.z() >= h_) return DBL_MAX;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if (y <= -w_ / 2 || y >= w_ / 2) return DBL_MAX;
				else if (y >= fsy2_ || y <= fsy1_) return dist;
			}

			// Транспорт внутри закутка
			double d = getDistanceToMlcRectangleConeInside(pp, u, -l_, -l_ + fsx1_, fsy1_, fsy2_);
			if (d == DBL_MAX) return DBL_MAX;
			// В закутке мым могли просто стартовать. 
			// Поэтому остается вариант вылета из него через открытую боковую поверхность.
			pp += u * d;
			if (pp.x() < -l_) return DBL_MAX;
			else return d + dist;
		}
	}
	else
	{
		// Сначала переводим за плоскость выступающей части
		if (fsx2_ > 0)	// выступает центральная часть 
		{
			if (x >= l_ + fsx2_)
			{
				if (u.x() >= 0) return DBL_MAX;
				double d = (fsx2_ + l_ - x) / u.x();
				pp += u * d;
				dist += d;
				if (pp.z() <= 0 || pp.z() >= h_) return DBL_MAX;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if ((y <= -w_ / 2 && u.y() <= 0) || (y >= w_ / 2 && u.y() >= 0)) return DBL_MAX;
				else if (y >= fsy1_ && y <= fsy2_) return dist;
			}
			// Транспорт в двух угловых закутках
			if (y < 0)
			{
				double d = getDistanceToMlcRectangleConeInside(pp, u, l_, l_ + fsx2_, -w_ / 2, fsy1_);
				if (d == DBL_MAX) return DBL_MAX;
				pp += u * d;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if (y <= -w_ / 2) return DBL_MAX;
				else return d + dist;
			}
			else
			{
				double d = getDistanceToMlcRectangleConeInside(pp, u, l_, l_ + fsx2_, fsy2_, w_ / 2);
				if (d == DBL_MAX) return DBL_MAX;
				pp += u * d;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if (y >= w_ / 2) return DBL_MAX;
				else return d + dist;
			}
		}
		else
		{
			if (x >= l_)
			{
				if (u.x() <= 0) return DBL_MAX;
				double d = (-l_ - x) / u.x();
				pp += u * d;
				dist += d;
				if (pp.z() <= 0 || pp.z() >= h_) return DBL_MAX;
				y = pp.y() * focus_ / (focus_ - pp.z());
				if (y <= -w_ / 2 || y >= w_ / 2) return DBL_MAX;
				else if (y >= fsy2_ || y <= fsy1_) return dist;
			}

			// Транспорт внутри закутка
			double d = getDistanceToMlcRectangleConeInside(pp, u, l_, l_ + fsx2_, fsy1_, fsy2_);
			if (d == DBL_MAX) return DBL_MAX;
			pp += u * d;
			if (pp.x() >= l_) return DBL_MAX;
			else return d + dist;
		}
	}
	return DBL_MAX;
}

double mcTransportMLC::getDistanceToMlcInsideField(const geomVector3D& p, const geomVector3D& u) const
{
	// Перебираем выходы из блока пока не вылетим, не столкнемся с границей или не перейдем в другой блок
	double dist = 0;
	double y = p.y() * focus_ / (focus_ - p.z());
	geomVector3D pp(p);

	// Цикл нужен для перебора пока не вывалимся.
	// На всякий случай не бесконечный.
	for (int k = 0; k < 3; k++)
	{
		if (y <= fsy1_)
		{
			double d = getDistanceToMlcInsideFieldBlock(pp, u, 0, 0, -w_ / 2, fsy1_);
			if (d == DBL_MAX) return DBL_MAX;
			pp += u * d;
			dist += d;
			y = pp.y() * focus_ / (focus_ - pp.z());
			if (fabs(y - fsy1_) <= DBL_EPSILON)
			{
				y = fsy1_;
				// Проверяем не уперлись ли в лепесток сбоку
				if (isOnOnTheLeftLeaf(pp, y, fsx1_) || isOnOnTheLeftLeaf(geomVector3D(-pp.x(), pp.y(), pp.z()), y, -fsx2_))
					break;
			}
			else if (fabs(y + w_/2) <= DBL_EPSILON) return DBL_MAX;
			else break;		// попали в торец лепестков
		}
		if (y >= fsy1_ && y < fsy2_)
		{
			double d = getDistanceToMlcInsideFieldBlock(pp, u, fsx1_, fsx2_, fsy1_, fsy2_);
			if (d == DBL_MAX) return DBL_MAX;
			pp += u * d;
			dist += d;
			y = pp.y() * focus_ / (focus_ - pp.z());
			if (fabs(y - fsy2_) <= DBL_EPSILON)
			{
				y = fsy2_;
				if (isOnOnTheLeftLeaf(pp, y, fsx1_) || isOnOnTheLeftLeaf(geomVector3D(-pp.x(), pp.y(), pp.z()), y, 0))
					break;
			}
			else if (fabs(y - fsy1_) <= DBL_EPSILON)
			{
				y = fsy1_;
				if (isOnOnTheLeftLeaf(pp, y, fsx1_) || isOnOnTheLeftLeaf(geomVector3D(-pp.x(), pp.y(), pp.z()), y, 0))
					break;
			}
			else break;
		}
		if (y >= fsy2_)
		{
			double d = getDistanceToMlcInsideFieldBlock(pp, u, 0, 0, fsy2_, w_ / 2);
			if (d == DBL_MAX) return DBL_MAX;
			pp += u * d;
			dist += d;
			y = pp.y() * focus_ / (focus_ - pp.z());
			if (fabs(y - fsy2_) <= DBL_EPSILON)
			{
				y = fsy2_;
				if (isOnOnTheLeftLeaf(pp, y, fsx1_) || isOnOnTheLeftLeaf(geomVector3D(-pp.x(), pp.y(), pp.z()), y, -fsx2_))
					break;
			}
			else if (fabs(y - w_/2) <= DBL_EPSILON) return DBL_MAX;
			else break;		// попали в торец лепестков
		}
	}

	return dist;
}

double mcTransportMLC::getDistanceToMlcInsideFieldBlock(const geomVector3D& p, const geomVector3D& u, double x1, double x2, double y1, double y2) const
{
	double vx = u.x(), vy = u.y(), vz = u.z();
	double x = p.x(), y = p.y(), z = p.z();

	// Из-за вогнутости поверхности нужно сначада проверять столкновение со скругленными торцами
	geomVector3D vv(vx, vz, vy);
	geomVector3D pp(x - x1 + r_, z - h_ / 2, 0);
	double d1 = mcGeometry::getDistanceToInfiniteCylinderOutside(pp, vv, r_);

	pp.set(x - r_ - x2, z - h_ / 2, 0);
	double d2 = mcGeometry::getDistanceToInfiniteCylinderOutside(pp, vv, r_);

	double dist = MIN(d1, d2);

	// Места пересечения боковых поверхностей
	double cd1 = DBL_MAX, cd2 = DBL_MAX;

	geomVector3D zxn1(0, focus_, y1), zxn2(0, -focus_, -y2);
	zxn1.normalize();
	zxn2.normalize();
	geomVector3D py1(0, y1, 0), py2(0, y2, 0);

	double dir = -(u * zxn1);
	if (dir > 0)
		cd1 = ((p - py1) * zxn1) / dir;

	dir = -(u * zxn2);
	if (dir > 0)
		cd2 = ((p - py2) * zxn2) / dir;

	dist = MIN(dist, MIN(cd1, cd2));
	if (dist == DBL_MAX) return DBL_MAX;

	// Проверяем, произошло ли пересечение в слое MLC
	double dd = z + vz * dist;
	if (dd <= 0 || dd >= h_) return DBL_MAX;
	else return dist;
}

double mcTransportMLC::getDistanceToMlcRectangleConeInside(const geomVector3D& p, const geomVector3D& u, double x1, double x2, double y1, double y2) const
{
	double x = p.x();
	double ux = u.x();
	double cd1 = DBL_MAX, cd2 = DBL_MAX, cd3 = DBL_MAX;

	if (ux > 0)
		cd1 = (x2 - x) / ux;
	else if (ux < 0)
		cd1 = (x1 - x) / ux;

	geomVector3D zxn1(0, focus_, y1), zxn2(0, -focus_, -y2);
	zxn1.normalize();
	zxn2.normalize();
	geomVector3D py1(0, y1, 0), py2(0, y2, 0);

	double dir = -(u * zxn1);
	if (dir > 0)
		cd2 = ((p - py1) * zxn1) / dir;

	dir = -(u * zxn2);
	if (dir > 0)
		cd3 = ((p - py2) * zxn2) / dir;

	double dist = MIN(cd1, MIN(cd2, cd3));
	if (dist == DBL_MAX) return DBL_MAX;
	double z = p.z() + dist * u.z();
	return (z >= 0 && z <= h_) ? dist : DBL_MAX;
}

bool mcTransportMLC::isInZPlaneHit(double x, double y) const
{
	if ((y >= -w_ / 2 && y <= fsy1_) || (y <= w_ / 2 && y >= fsy2_))
	{
		if (x < -l_ || x > l_) return false;
		else if (x > -dx_ && x < dx_) return false;
		else return true;
	}
	else if (y >= fsy1_ && y <= fsy2_)
	{
		if (x < (-l_ - fsx1_) || x >(l_ + fsx2_)) return false;
		else if (x > (-dx_ + fsx1_) && x < (dx_ + fsx2_)) return false;
		else return true;
	}
	return false;
}

bool mcTransportMLC::isOnOutsideLeaf(const geomVector3D& p) const
{
	double x = p.x();
	double z = p.z();

	if (x < -l_ || x > l_) return false;
	else if (x <= 0 && x > -r_)
	{
		double xx = r_ + x;
		if ((xx * xx + z * z) > r2_)
			return false;
	}
	else if (x >= 0 && x < r_)
	{
		double xx = r_ - x;
		if ((xx * xx + z * z) > r2_)
			return false;
	}

	return true;
}

bool mcTransportMLC::isOnOnTheLeftLeaf(const geomVector3D& p, double y0, double xx) const
{
	double z = p.z();
	if (z < 0 || z > h_) return false;

	double y = p.y() * focus_ / (focus_ - z);
	if (fabs(y - y0) > DBL_EPSILON) return false;

	double x = p.x();
	if(x < xx - l_) return false;

	double d = x + r_ - xx;
	if (d > 0 && (d * d + z * z) > r2_)
		return false;
	else
		return true;
}

bool mcTransportMLC::isOutsideMLC(double x, double y) const
{
	// Требуется только определеить одно из возможных положений -
	// внутри открытой части коллиматора или за пределами коллиматор
	if (y <= -w_ / 2 || y >= w_ / 2 || x <= -l_ || x >= l_)
		return true;
	else
		return false;
}

void mcTransportMLC::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Focus:\t" << focus_ << endl;
	os << "R:\t" << r_ << endl;
	os << "H:\t" << h_ << endl;
	os << "Width:\t" << w_ << endl;
	os << "Length:\t" << l_ << endl;
	os << "FSX1:\t" << fsx1_ << endl;
	os << "FSX2:\t" << fsx2_ << endl;
	os << "FSY1:\t" << fsy1_ << endl;
	os << "FSY2:\t" << fsy2_ << endl;
}

void mcTransportMLC::dumpVRML(ostream& os) const
{
	os << "# Polygon side hole: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	dumpVRMLMLC(os);

	os << "  ]" << endl;
	os << "}" << endl;
}

void mcTransportMLC::dumpVRMLMLC(ostream& os) const
{
	int i = 0, k, np = 24 * 2 + 6 * 7 * 2;
	vector<geomVector3D> p(np);

	double xx[8], yy[4], ss[7];
	double dh = h_ / 6;
	double dx = r_ * (1 - sqrt(1.0 - 9.0*dh*dh / r2_));
	double sf = (focus_ - h_) / focus_;		// масштаб перевода координаты y на верхнюю поверхность

	xx[0] = fsx1_ - l_;
	xx[1] = -l_;
	xx[2] = fsx1_ - dx;
	xx[3] = -dx;
	xx[4] = dx;
	xx[5] = fsx2_ + dx;
	xx[6] = l_;
	xx[7] = fsx2_ + l_;

	yy[0] = -w_ / 2;
	yy[1] = fsy1_;
	yy[2] = fsy2_;
	yy[3] = w_ / 2;

	for (k = 0; k < 7; k++)
		ss[k] = dx - r_ * (1 - sqrt(1.0 - (k - 3) * (k - 3) * dh * dh / r2_));

	// Верхняя и нижняя поверхности коллиматора
	p[i++] = geomVector3D(xx[1], yy[3], 0) * mttow_;
	p[i++] = geomVector3D(xx[1], yy[3] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[3], 0) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[3] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[3], 0) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[3] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[3], 0) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[3] * sf, h_) * mttow_;

	p[i++] = geomVector3D(xx[0], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[0], yy[2] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[1], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[1], yy[2] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[2], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[2], yy[2] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[2] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[2] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[5], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[5], yy[2] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[2] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[7], yy[2], 0) * mttow_;
	p[i++] = geomVector3D(xx[7], yy[2] * sf, h_) * mttow_;

	p[i++] = geomVector3D(xx[0], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[0], yy[1] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[1], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[1], yy[1] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[2], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[2], yy[1] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[1] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[1] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[5], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[5], yy[1] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[1] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[7], yy[1], 0) * mttow_;
	p[i++] = geomVector3D(xx[7], yy[1] * sf, h_) * mttow_;

	p[i++] = geomVector3D(xx[1], yy[0], 0) * mttow_;
	p[i++] = geomVector3D(xx[1], yy[0] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[0], 0) * mttow_;
	p[i++] = geomVector3D(xx[3], yy[0] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[0], 0) * mttow_;
	p[i++] = geomVector3D(xx[4], yy[0] * sf, h_) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[0], 0) * mttow_;
	p[i++] = geomVector3D(xx[6], yy[0] * sf, h_) * mttow_;

	int lcount = i;		// счетчик точек, отвечающих за верхнюю и нижнюю поверхности

	// Точки скругленных поверхностей МЛК
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[3] + ss[k], yy[3] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[3] + ss[k], yy[2] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[2] + ss[k], yy[2] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[2] + ss[k], yy[1] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[3] + ss[k], yy[1] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[3] + ss[k], yy[0] * (1 - dh * k / focus_), dh * k) * mttow_;

	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[4] - ss[k], yy[3] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[4] - ss[k], yy[2] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[5] - ss[k], yy[2] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[5] - ss[k], yy[1] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[4] - ss[k], yy[1] * (1 - dh * k / focus_), dh * k) * mttow_;
	for (k = 0; k < 7; k++) p[i++] = geomVector3D(xx[4] - ss[k], yy[0] * (1 - dh * k / focus_), dh * k) * mttow_;

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

	// Выводим все координаты
	for (i = 0; i < np; i++) {
		os << "                    " << p[i].x() << ' ' << p[i].y() << ' ' << p[i].z();
		if (i < np - 1) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	// Нижняя и верхняя поверхности
	//os << "                0, 10, 8, 24, 26, 40, 42, 30, 28, 12, 14, 2, -1," << endl;
	//os << "                1, 3, 15, 13, 29, 31, 43, 41, 27, 25, 9, 11, -1," << endl;

	//os << "                4, 6, 20, 22, 38, 36, 46, 44, 32, 34, 18, 16, -1," << endl;
	//os << "                5, 17, 19, 35, 33, 45, 47, 37, 39, 23, 21, 7, -1," << endl;

	os << "                0, 2, 14, 10, -1," << endl;
	os << "                1, 11, 15, 3, -1," << endl;
	os << "                4, 6, 20, 16, -1," << endl;
	os << "                5, 17, 21, 7, -1," << endl;

	os << "                8, 12, 28, 24, -1," << endl;
	os << "                9, 25, 29, 13, -1," << endl;
	os << "                18, 22, 38, 34, -1," << endl;
	os << "                19, 35, 39, 23, -1," << endl;

	os << "                26, 30, 42, 40, -1," << endl;
	os << "                27, 41, 43, 31, -1," << endl;
	os << "                32, 36, 46, 44, -1," << endl;
	os << "                33, 45, 47, 37, -1," << endl;

	// Внешние торцы
	os << "                0, 10, 11, 1, -1," << endl;
	os << "                10, 8, 9, 11, -1," << endl;
	os << "                8, 24, 25, 9, -1," << endl;
	os << "                24, 26, 27, 25, -1," << endl;
	os << "                26, 40, 41, 27, -1," << endl;

	os << "                20, 6, 7, 21, -1," << endl;
	os << "                22, 20, 21, 23, -1," << endl;
	os << "                38, 22, 23, 39, -1," << endl;
	os << "                36, 38, 39, 37, -1," << endl;
	os << "                46, 36, 37, 47, -1," << endl;

	// Скругленные боковые поверхности
	os << "                0, 1, ";
	for (k = 6; k >= 0; k--) os << (lcount + k) << ", ";
	os << " -1," << endl;

	os << "                7, 6, ";
	for (k = 0; k < 7; k++) os << (lcount + 6 * 7 + k) << ", ";
	os << " -1," << endl;

	os << "                41, 40, ";
	for (k = 0; k < 7; k++) os << (lcount + 5 * 7 + k) << ", ";
	os << " -1," << endl;

	os << "                46, 47, ";
	for (k = 6; k >= 0; k--) os << (lcount + 11 * 7 + k) << ", ";
	os << " -1," << endl;

	// Внутренние боковые поверхности
	os << "                ";
	for (k = 0; k < 7; k++) os << (lcount + 1 * 7 + k) << ", ";
	for (k = 6; k >= 0; k--) os << (lcount + 2 * 7 + k) << ", ";
	os << " -1," << endl;

	os << "                ";
	for (k = 0; k < 7; k++) os << (lcount + 8 * 7 + k) << ", ";
	for (k = 6; k >= 0; k--) os << (lcount + 7 * 7 + k) << ", ";
	os << " -1," << endl;

	os << "                ";
	for (k = 0; k < 7; k++) os << (lcount + 3 * 7 + k) << ", ";
	for (k = 6; k >= 0; k--) os << (lcount + 4 * 7 + k) << ", ";
	os << " -1," << endl;

	os << "                ";
	for (k = 0; k < 7; k++) os << (lcount + 10 * 7 + k) << ", ";
	for (k = 6; k >= 0; k--) os << (lcount + 9 * 7 + k) << ", ";
	os << " -1," << endl;

	// Круглые торцы лепестков
	os << "                ";
	for (k = 0; k < 6; k++) os << (lcount + 1 * 7 + k) << ", " << (lcount + 0 * 7 + k) << ", " << (lcount + 0 * 7 + k + 1) << ", " << (lcount + 1 * 7 + k + 1) << " -1," << endl;
	os << "                ";
	for (k = 0; k < 6; k++) os << (lcount + 3 * 7 + k) << ", " << (lcount + 2 * 7 + k) << ", " << (lcount + 2 * 7 + k + 1) << ", " << (lcount + 3 * 7 + k + 1) << " -1," << endl;
	os << "                ";
	for (k = 0; k < 6; k++) os << (lcount + 5 * 7 + k) << ", " << (lcount + 4 * 7 + k) << ", " << (lcount + 4 * 7 + k + 1) << ", " << (lcount + 5 * 7 + k + 1) << " -1," << endl;

	os << "                ";
	for (k = 0; k < 6; k++) os << (lcount + 6 * 7 + k) << ", " << (lcount + 7 * 7 + k) << ", " << (lcount + 7 * 7 + k + 1) << ", " << (lcount + 6 * 7 + k + 1) << " -1," << endl;
	os << "                ";
	for (k = 0; k < 6; k++) os << (lcount + 8 * 7 + k) << ", " << (lcount + 9 * 7 + k) << ", " << (lcount + 9 * 7 + k + 1) << ", " << (lcount + 8 * 7 + k + 1) << " -1," << endl;
	os << "                ";
	for (k = 0; k < 6; k++) os << (lcount + 10 * 7 + k) << ", " << (lcount + 11 * 7 + k) << ", " << (lcount + 11 * 7 + k + 1) << ", " << (lcount + 10 * 7 + k + 1) << " -1," << endl;

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}
