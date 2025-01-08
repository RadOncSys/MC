#include "mcTransportMantleBlock.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportMantleBlock::mcTransportMantleBlock(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, 
	double r1, double h, std::vector<double> plgnX, std::vector<double> plgnY)
	: mcTransport(orgn, z, x), x_(plgnX), y_(plgnY), r1_(r1), h_(h)
{
	nsides_ = plgnX.size();
	geomVector3D vy(0, 0, 1);
	sides_ = std::make_unique<std::vector<mcGeomRectSide>>(nsides_);

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
	if (pz <= 0)
	{
		dist = -pz / uz;
		c += p.u * dist;
		if (c.lengthXY() > r1_)
			return DBL_MAX;
	}
	else if (pz >= h_)
	{
		dist = -(pz - h_) / uz;
		c += p.u * dist;
		if (c.lengthXY() > r1_)
			return DBL_MAX;
	}

	// Если мы переместились на поверхность, то отсекаем ситуацию что мы уже в теле объекта

	if (dist > 0)
	{
		if (!isPointInPlgn(c.x(), c.y()))
			return dist;
		else
		{
			double cd2 = DBL_MAX;
			if (pz >= 0 && pz <= h_)
			{
				for (int i = 0; i < nsides_; i++)
				{
					double d = sides_->at(i).getDistance(p.p, p.u);
					if (d < cd2) cd2 = d;
				}
			}
			return cd2;
		}
	}
}

bool mcTransportMantleBlock::isPointInPlgn(double x, double y) const
{
	// Собираем все сечения полигона линиями y = const, где y - координата тестируемой точки.

	// Счетчик пересечений слева от точки
	int count = 0;

	for (int i = 0; i < nsides_; i++)
	{
		int i1 = (i + 1) % nsides_;
		double dy = y_[i] - y;
		double dy1 = y_[i1] - y;
		if (dy * dy1 < 0)
		{
			double cx = x_[i] + (x_[i1] - x_[i]) * (y - y_[i]) / (y_[i1] - y_[i]);
			if (cx < x)
				count++;
		}

		// TODO: Есть проблема, когда y совпадает с координатой одной или более точек полигона.


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
	os << "# Conical hole: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	dumpVRMLCylinder(os, r1_, 0, h_, 0, 0);

	os << "  ]" << endl;
	os << "}" << endl;
}
