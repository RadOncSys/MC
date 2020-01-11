#include "mcTransportRectanglePolygonSideHole.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportRectanglePolygonSideHole::mcTransportRectanglePolygonSideHole(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
	double dx, double dy, const std::vector<double>& z, const std::vector<double>& x, const std::vector<double>& y)
	:mcTransport(orgn, vz, vx), dx_(dx), dy_(dy), z_(z), x_(x), y_(y), fsx1_(0), fsx2_(0), fsy1_(0), fsy2_(0)
{
	nlayers_ = (int)z_.size() - 1;
	xmax_ = x_[0];
	ymax_ = y_[0];
	zmax_ = z_.back();

	// Наклоны граней 
	cosx_.resize(nlayers_);
	sinx_.resize(nlayers_);
	cosy_.resize(nlayers_);
	siny_.resize(nlayers_);

	for (int i = 0; i < nlayers_; i++)
	{
		double a = x_[i + 1] - x_[i], b = z_[i + 1] - z_[i], c = sqrt(a * a + b * b);
		sinx_[i] = a / c;
		cosx_[i] = b / c;
		a = y_[i + 1] - y_[i], c = sqrt(a * a + b * b);
		siny_[i] = a / c;
		cosy_[i] = b / c;
		xmax_ = MAX(xmax_, x_[i + 1]);
		ymax_ = MAX(ymax_, y_[i + 1]);
	}
}

mcTransportRectanglePolygonSideHole::~mcTransportRectanglePolygonSideHole(void)
{
}

double mcTransportRectanglePolygonSideHole::getDistanceInside(mcParticle& p) const
{
	double z = p.p.z(), vz = p.u.z();

	// Проверяем пересечение с внутренней конусообразной поверхностью имеет смысл рассматривать в первую очередь,
	// так как она внутри и, если пересечение существует, то оно ближайшее
	for (int i = 0; i < nlayers_; i++)
	{
		double cd = mcGeometry::getDistanceToRectangleConeOutside(p.p, p.u,
			-x_[i] + fsx1_, x_[i] + fsx2_, cosx_[i], sinx_[i],
			-y_[i] + fsy1_, y_[i] + fsy2_, cosy_[i], siny_[i], z_[i]);

		if (cd == DBL_MAX) continue;
		geomVector3D c = p.p + (p.u * cd);
		if (c.z() >= z_[i] && c.z() <= z_[i + 1])	// Пересечение внутри слоя
			return cd;
	}

	// Проверяем прохождение через переднюю или заднюю поверхность
	// (интерференция с боковыми стенками отсутствует, так как без внутренней дырки,
	// рассмотренной выше, внутренний объем является выпуклым, т.е. 
	// между любыми двумя точками принадлежащими объему, нет пересечений)
	if (vz < 0)
	{
		double cd = -z / vz;
		geomVector3D c = p.p + (p.u * cd);
		double x = c.x(), y = c.y();
		if (x >= (-dx_ + fsx1_) && x <= (dx_ + fsx2_) && y >= (-dy_ + fsy1_) && y <= (dy_ + fsy2_) &&
			(x <= (-x_.front() + fsx1_) || x >= (x_.front() + fsx2_) || y <= (-y_.front() + fsy1_) || y >= (y_.front() + fsy2_)))
			return cd;
	}
	else if (vz > 0)
	{
		double cd = (zmax_ - z) / vz;
		geomVector3D c = p.p + (p.u * cd);
		double x = c.x(), y = c.y();
		if (x >= (-dx_ + fsx1_) && x <= (dx_ + fsx2_) && y >= (-dy_ + fsy1_) && y <= (dy_ + fsy2_) &&
			(x <= (-x_.back() + fsx1_) || x >= (x_.back() + fsx2_) || y <= (-y_.back() + fsy1_) || y >= (y_.back() + fsy2_)))
			return cd;
	}

	// Проверяем пересечение с внешними боковыми гранями
	double cd = mcGeometry::getDistanceToRectanglePipeInside(p.p, p.u, -dx_ + fsx1_, dx_ + fsx2_, -dy_ + fsy1_, dy_ + fsy2_);
	if (cd != DBL_MAX)
	{
		geomVector3D c = p.p + (p.u * cd);
		if (c.z() >= 0 && c.z() <= zmax_)
			return cd;
	}

	// Important!
	// Хороший прием поиска проблемы - повторение расета для той же частицы.

	// Отладка: Повторить расчет, чтобы понять причину.
	// Не допускать в норме, так как приведет к бесконечной рекурсии
	//return getDistanceInside(p.p, p.u);

	//throw std::exception("mcTransportRectanglePolygonSideHole::getDistanceInside: particle inside should have distance");
	return DBL_MAX;
}

double mcTransportRectanglePolygonSideHole::getDistanceOutside(mcParticle& p) const
{
	double z = p.p.z(), vz = p.u.z();
	geomVector3D c(p.p);	// Требуется для проекции на переднюю или заднюю область объекта

	// Частица перед объектом
	if (z <= DBL_EPSILON)
	{
		if (vz <= 0) return DBL_MAX;

		// Доводим частицу до ближайшей плоскости объекта
		double cd = fabs(z / vz); // fabs - борьба с погрешностями расчетов, когда частица на поверхности, но не точно
		c = p.p + (p.u * cd);
		double x = c.x(), y = c.y();

		// Если частица окажется на поверхности объекта, то это уже пересечение
		if (x >= -dx_ + fsx1_ && x <= dx_ + fsx2_ && y >= -dy_ + fsy1_ && y <= dy_ + fsy2_ &&
			!(x > -x_.front() + fsx1_ && x < x_.front() + fsx2_ &&
				y > -y_.front() + fsx1_ && y < y_.front() + fsy2_))
			return cd;
		else
		{
			double sideDistance = getDistanceOutsideWithinLayer(c, p.u);
			return sideDistance == DBL_MAX ? sideDistance : cd + sideDistance;
		}
	}

	// Частица за объектом
	else if (z >= zmax_ - DBL_EPSILON)
	{
		if (vz >= 0) return DBL_MAX;
		double cd = fabs((zmax_ - z) / vz);
		c = p.p + (p.u * cd);
		double x = c.x(), y = c.y();

		// Если частица окажется на поверхности объекта, то это уже пересечение
		if (x >= -dx_ + fsx1_ && x <= dx_ + fsx2_ && y >= -dy_ + fsy1_ && y <= dy_ + fsy2_ &&
			!(x > -x_.back() + fsx1_ && x < x_.back() + fsx2_ &&
				y > -y_.back() + fsx1_ && y < y_.back() + fsy2_))
			return cd;
		else
		{
			double sideDistance = getDistanceOutsideWithinLayer(c, p.u);
			return sideDistance == DBL_MAX ? sideDistance : cd + sideDistance;
		}
	}

	else
		return getDistanceOutsideWithinLayer(c, p.u);
}

double mcTransportRectanglePolygonSideHole::getDistanceOutsideWithinLayer(const geomVector3D& p, const geomVector3D& v) const
{
	// Если частица за пределами периметра, то пересечение возможно только с внешней вертикальной гранью.
	// В дырке тащим ее по слоям.
	double x = p.x(), y = p.y();

	if (x <= -dx_ + fsx1_ || x >= dx_ + fsx2_ || y <= -dy_ + fsy1_ || y >= dy_ + fsy2_)
	{
		double cd = mcGeometry::getDistanceToRectanglePipeOutside(p, v, -dx_ + fsx1_, dx_ + fsx2_, -dy_ + fsy1_, dy_ + fsy2_);
		if (cd == DBL_MAX)
			return DBL_MAX;

		// Точка пересечения с прямоугольным цилиндром, определяющим внешние грани
		auto c = p + (v * cd);
		if (c.z() <= 0 || c.z() >= zmax_)
			return DBL_MAX;
		else
			return cd;
	}
	else
	{
		// Оставшиеся частицы могут быть только в дырке. Проверяем все грани.
		for (int i = 0; i < nlayers_; i++)
		{
			double cd = mcGeometry::getDistanceToRectangleConeInside(p, v,
				-x_[i] + fsx1_, x_[i] + fsx2_, cosx_[i], sinx_[i],
				-y_[i] + fsy1_, y_[i] + fsy2_, cosy_[i], siny_[i], z_[i]);
			if (cd == DBL_MAX) continue;
			auto c = p + (v * cd);
			if (c.z() >= z_[i] && c.z() <= z_[i + 1])
				return cd;	// Пересечение внутри слоя
		}
		return DBL_MAX;
	}
}

double mcTransportRectanglePolygonSideHole::getDNearInside(const geomVector3D& p) const
{
	// Вычисление расстояние нужно только для ускорения транспорта частиц.
	// Если всегда возвращать 0, то не будет ускорения за счет сокращения обращений к 
	// функции вычисления пересечений с границами.
	// В данном случае вычисление минимального расстояния затруднительно.
	// Поэтому берем только максимальное вписанное прямоугольное кольцо.
	// В остальных частях считаем расстояние до границы нулевым.

	double x = p.x(), y = p.y(), z = p.z();
	if (x <= (xmax_ + fsx2_) && x >= (-xmax_ + fsx1_) && y <= (ymax_ + fsy2_) && y >= (-ymax_ + fsy1_))
		return 0;

	double dd = MIN(MIN(ZORP(x + dx_ - fsx1_), ZORP(dx_ + fsx2_ - x)), MIN(ZORP(y + dy_ - fsy1_), ZORP(dy_ + fsy2_ - y)));
	dd = MIN(dd, ZORP(z));
	dd = MIN(dd, ZORP(zmax_ - z));

	double dx1 = (fsx1_ - xmax_) - x, dx2 = x - (fsx2_ + xmax_);
	double dy1 = (fsy1_ - ymax_) - x, dy2 = x - (fsy2_ + ymax_);


	if (dx1 > 0)
	{
		if (dy1 > 0)
			return MIN(dd, sqrt(dx1 * dx1 + dy1 * dy1));
		else if (dy2 > 0)
			return MIN(dd, sqrt(dx1 * dx1 + dy2 * dy2));
		else
			return MIN(dd, dx1);
	}
	else if (dx2 > 0)
	{
		if (dy1 > 0)
			return MIN(dd, sqrt(dx2 * dx2 + dy1 * dy1));
		else if (dy2 > 0)
			return MIN(dd, sqrt(dx2 * dx2 + dy2 * dy2));
		else
			return MIN(dd, dx2);
	}
	else
	{
		if (dy1 > 0)
			return MIN(dd, dy1);
		else
			return MIN(dd, dy2);
	}
}

void mcTransportRectanglePolygonSideHole::SetFieldSize(double x1, double x2, double y1, double y2)
{
	fsx1_ = x1;
	fsx2_ = x2;
	fsy1_ = y1;
	fsy2_ = y2;
}

void mcTransportRectanglePolygonSideHole::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "DX:\t" << dx_ << endl;
	os << "DY:\t" << dy_ << endl;
	os << "NLayers:\t" << nlayers_ << endl;
	os << "XMax:\t" << xmax_ << endl;
	os << "YMax:\t" << ymax_ << endl;
	os << "ZMax:\t" << zmax_ << endl;
	os << "FsX1:\t" << fsx1_ << endl;
	os << "FsX2:\t" << fsx2_ << endl;
	os << "FsY1:\t" << fsy1_ << endl;
	os << "FsY2:\t" << fsy2_ << endl;

	unsigned i;
	os << "Polygon:" << endl;
	os << "X:";
	for (i = 0; i < x_.size(); i++) os << "\t" << x_[i];
	os << endl;

	os << "Y:";
	for (i = 0; i < y_.size(); i++) os << "\t" << y_[i];
	os << endl;

	os << "Z:";
	for (i = 0; i < z_.size(); i++) os << "\t" << z_[i];
	os << endl;

	os << "Cosines:" << endl;
	os << "SinX:";
	for (i = 0; i < sinx_.size(); i++) os << "\t" << sinx_[i];
	os << endl;
	os << "CosX:";
	for (i = 0; i < cosx_.size(); i++) os << "\t" << cosx_[i];
	os << endl;
	os << "SinY:";
	for (i = 0; i < siny_.size(); i++) os << "\t" << siny_[i];
	os << endl;
	os << "CosY:";
	for (i = 0; i < cosy_.size(); i++) os << "\t" << cosy_[i];
	os << endl;
}

void mcTransportRectanglePolygonSideHole::dumpVRML(ostream& os) const
{
	os << "# Polygon side hole: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	dumpVRMLPolygonSideHole(os);

	os << "  ]" << endl;
	os << "}" << endl;
}

void mcTransportRectanglePolygonSideHole::dumpVRMLPolygonSideHole(ostream& os) const
{
	unsigned i = 0, np = 16 + (unsigned)z_.size() * 4, k;
	vector<geomVector3D> p(np);

	// Передняя поверхность
	p[i++] = geomVector3D(-dx_ + fsx1_, -dy_ + fsy1_, 0) * mttow_;
	p[i++] = geomVector3D(-dx_ + fsx1_, dy_ + fsy2_, 0) * mttow_;
	p[i++] = geomVector3D(dx_ + fsx2_, dy_ + fsy2_, 0) * mttow_;
	p[i++] = geomVector3D(dx_ + fsx2_, -dy_ + fsy1_, 0) * mttow_;
	p[i++] = geomVector3D(-x_[0] + fsx1_, -y_[0] + fsy1_, 0) * mttow_;
	p[i++] = geomVector3D(-x_[0] + fsx1_, y_[0] + fsy2_, 0) * mttow_;
	p[i++] = geomVector3D(x_[0] + fsx2_, y_[0] + fsy2_, 0) * mttow_;
	p[i++] = geomVector3D(x_[0] + fsx2_, -y_[0] + fsy1_, 0) * mttow_;

	// Задняя поверхность
	p[i++] = geomVector3D(-dx_ + fsx1_, -dy_ + fsy1_, zmax_) * mttow_;
	p[i++] = geomVector3D(-dx_ + fsx1_, dy_ + fsy2_, zmax_) * mttow_;
	p[i++] = geomVector3D(dx_ + fsx2_, dy_ + fsy2_, zmax_) * mttow_;
	p[i++] = geomVector3D(dx_ + fsx2_, -dy_ + fsy1_, zmax_) * mttow_;
	p[i++] = geomVector3D(-x_.back() + fsx1_, -y_.back() + fsy1_, zmax_) * mttow_;
	p[i++] = geomVector3D(-x_.back() + fsx1_, y_.back() + fsy2_, zmax_) * mttow_;
	p[i++] = geomVector3D(x_.back() + fsx2_, y_.back() + fsy2_, zmax_) * mttow_;
	p[i++] = geomVector3D(x_.back() + fsx2_, -y_.back() + fsy1_, zmax_) * mttow_;

	// Дырка
	for (k = 0; k < z_.size(); k++)
	{
		p[i++] = geomVector3D(-x_[k] + fsx1_, -y_[k] + fsy1_, z_[k]) * mttow_;
		p[i++] = geomVector3D(-x_[k] + fsx1_, y_[k] + fsy2_, z_[k]) * mttow_;
		p[i++] = geomVector3D(x_[k] + fsx2_, y_[k] + fsy2_, z_[k]) * mttow_;
		p[i++] = geomVector3D(x_[k] + fsx2_, -y_[k] + fsy1_, z_[k]) * mttow_;
	}

	// Поверхность дырки (точки на основаниях повторяются для прозрачности)

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

	for (i = 0; i < p.size(); i++) {
		os << "                    " << p[i].x() << ' ' << p[i].y() << ' ' << p[i].z();
		if (i < 15) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	// Внешняя боковая поверхность
	os << "                0, 8, 9, 1, -1," << endl;
	os << "                1, 9, 10, 2, -1," << endl;
	os << "                2, 10, 11, 3, -1," << endl;
	os << "                3, 11, 8, 0, -1," << endl;

	// Нижний торец
	os << "                0, 1, 5, 4, -1," << endl;
	os << "                1, 2, 6, 5, -1," << endl;
	os << "                2, 3, 7, 6, -1," << endl;
	os << "                3, 0, 4, 7, -1," << endl;

	// Поверхность дырки
	for (k = 0; k < (unsigned)nlayers_; k++)
	{
		unsigned s = 16 + k * 4;
		os << "                " << s << ", " << s + 1 << ", " << s + 5 << ", " << s + 4 << ", -1," << endl;
		os << "                " << s + 1 << ", " << s + 2 << ", " << s + 6 << ", " << s + 5 << ", -1," << endl;
		os << "                " << s + 2 << ", " << s + 3 << ", " << s + 7 << ", " << s + 6 << ", -1," << endl;
		os << "                " << s + 3 << ", " << s << ", " << s + 4 << ", " << s + 7 << ", -1," << endl;
	}

	// Верхний торец
	os << "                8, 12, 13, 9, -1," << endl;
	os << "                9, 13, 14, 10, -1," << endl;
	os << "                10, 14, 15, 11, -1," << endl;
	os << "                11, 15, 12, 8, -1" << endl;

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}
