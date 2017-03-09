#include "mcTransportRectanglePolygonSideHole.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportRectanglePolygonSideHole::mcTransportRectanglePolygonSideHole(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
	double dx, double dy, std::vector<double>& z, std::vector<double>& x, std::vector<double>& y)
	:mcTransport(orgn, vz, vx), dx_(dx), dy_(dy), z_(z), x_(x), y_(y)
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
		double a = x_[i + 1] - x_[i], b = z_[i + 1] - z_[i], c = sqrt(a*a + b*b);
		sinx_[i] = a / c;
		cosx_[i] = b / c;
		a = y_[i + 1] - y_[i], c = sqrt(a*a + b*b);
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
		double cd = mcGeometry::getDistanceToRectangleConeOutside(p.p, p.u, x_[i], cosx_[i], sinx_[i], y_[i], cosy_[i], siny_[i], z_[i]);
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
		double xx = fabs(c.x()), yy = fabs(c.y());
		if (xx <= dx_ && yy <= dy_ && (xx >= x_[0] || yy >= y_[0]))
			return cd;
	}
	else if (vz > 0)
	{
		double cd = (zmax_ - z) / vz;
		geomVector3D c = p.p + (p.u * cd);
		double xx = fabs(c.x()), yy = fabs(c.y());
		if (xx <= dx_ && yy <= dy_ && (xx >= x_.back() || yy >= y_.back()))
			return cd;
	}

	// Проверяем пересечение с внешними боковыми гранями
	double cd = mcGeometry::getDistanceToRectanglePipeInside(p.p, p.u, -dx_, dx_, -dy_, dy_);
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

		// Если частица окажется на поверхности объекта, то это уже пересечение
		double xx = fabs(c.x()), yy = fabs(c.y());
		if (xx <= dx_ && yy <= dy_ && (xx >= x_[0] || yy >= y_[0]))
			return cd;
	}

	// Частица за объектом
	else if (z >= zmax_ - DBL_EPSILON)
	{
		if (vz >= 0) return DBL_MAX;
		double cd = fabs((zmax_ - z) / vz);
		c = p.p + (p.u * cd);
		double xx = fabs(c.x()), yy = fabs(c.y());
		if (xx <= dx_ && yy <= dy_ && (xx >= x_.back() || yy >= y_.back()))
			return cd;
	}

	// Если частица дотянула до этого места, то она на уровне объекта.
	// Если частица за пределами периметра, то пересечение возможно только с внешней вертикальной гранью.
	// В дырке тащим ее по слоям.

	if (fabs(c.x()) >= dx_ || fabs(c.y()) >= dy_)
	{
		double cd = mcGeometry::getDistanceToRectanglePipeOutside(p.p, p.u, -dx_, dx_, -dy_, dy_);
		if (cd == DBL_MAX)
			return DBL_MAX;

		// Точка пересечения с прямоугольным цилиндром, определяющим внешние грани
		c = p.p + (p.u * cd);
		if (c.z() <= 0 || c.z() >= zmax_)
			return DBL_MAX;
		else
			return cd;
	}

	// Оставшиеся частицы могут быть только в дырке. Проверяем все грани.
	for (int i = 0; i < nlayers_; i++)
	{
		double cd = mcGeometry::getDistanceToRectangleConeInside(p.p, p.u, x_[i], cosx_[i], sinx_[i], y_[i], cosy_[i], siny_[i], z_[i]);
		if (cd == DBL_MAX) continue;
		c = p.p + (p.u * cd);
		if (c.z() >= z_[i] && c.z() <= z_[i + 1])
			return cd;	// Пересечение внутри слоя
	}

	return DBL_MAX;
}

double mcTransportRectanglePolygonSideHole::getDNearInside(const geomVector3D& p) const
{
	// Вычисление расстояние нужно только для ускорения транспорта частиц.
	// Если всегда возвращать 0, то не будет ускорения за счет сокращения обращений к 
	// функции вычисления пересечений с границами.
	// В данном случае вычисление минимального расстояния затруднительно.
	// Поэтому берем только максимальное вписанное прямоугольное кольцо.
	// В остальных частях считаем расстояние до границы нулевым.

	double x = fabs(p.x()), y = fabs(p.y()), z = p.z();
	if (x <= xmax_ && y <= ymax_)
		return 0;

	double dd = MIN(ZORP(dx_ - x), ZORP(dy_ - y));
	dd = MIN(dd, ZORP(z));
	dd = MIN(dd, ZORP(zmax_ - z));

	double dx = ZORP(x - xmax_);
	double dy = ZORP(y - ymax_);

	if (x < xmax_)
		return MIN(dd, dy);
	else if (y < ymax_)
		return MIN(dd, dx);
	else
		return MIN(dd, sqrt(dx*dx + dy*dy));
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
	p[i++] = geomVector3D(-dx_, -dy_, 0) * mttow_;
	p[i++] = geomVector3D(-dx_, dy_, 0) * mttow_;
	p[i++] = geomVector3D(dx_, dy_, 0) * mttow_;
	p[i++] = geomVector3D(dx_, -dy_, 0) * mttow_;
	p[i++] = geomVector3D(-x_[0], -y_[0], 0) * mttow_;
	p[i++] = geomVector3D(-x_[0], y_[0], 0) * mttow_;
	p[i++] = geomVector3D(x_[0], y_[0], 0) * mttow_;
	p[i++] = geomVector3D(x_[0], -y_[0], 0) * mttow_;

	// Задняя поверхность
	p[i++] = geomVector3D(-dx_, -dy_, zmax_) * mttow_;
	p[i++] = geomVector3D(-dx_, dy_, zmax_) * mttow_;
	p[i++] = geomVector3D(dx_, dy_, zmax_) * mttow_;
	p[i++] = geomVector3D(dx_, -dy_, zmax_) * mttow_;
	p[i++] = geomVector3D(-x_.back(), -y_.back(), zmax_) * mttow_;
	p[i++] = geomVector3D(-x_.back(), y_.back(), zmax_) * mttow_;
	p[i++] = geomVector3D(x_.back(), y_.back(), zmax_) * mttow_;
	p[i++] = geomVector3D(x_.back(), -y_.back(), zmax_) * mttow_;

	// Дырка
	for (k = 0; k < z_.size(); k++)
	{
		p[i++] = geomVector3D(-x_[k], -y_[k], z_[k]) * mttow_;
		p[i++] = geomVector3D(-x_[k], y_[k], z_[k]) * mttow_;
		p[i++] = geomVector3D(x_[k], y_[k], z_[k]) * mttow_;
		p[i++] = geomVector3D(x_[k], -y_[k], z_[k]) * mttow_;
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
