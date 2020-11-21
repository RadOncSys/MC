#include "mcTransportJawPairFocused.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportJawPairFocused::mcTransportJawPairFocused(void)
	:mcTransport()
{
}

mcTransportJawPairFocused::mcTransportJawPairFocused(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x,
	double sad, double scd, double h, double dx, double dy)
	: mcTransport(orgn, z, x)
	, sad_(sad), scd_(scd), h_(h), dx_(dx), dy_(dy)
{
}

mcTransportJawPairFocused::~mcTransportJawPairFocused(void)
{
}

void mcTransportJawPairFocused::setFS(double x1, double x2) 
{
	fsx1_ = x1; 
	fsx2_ = x2;
	fscenter_ = 0.5 * (x1 + x2); 

	//x2l_ = fsx1_ * scd_ / sad_;
	x2l_ = fsx1_;
	x1l_ = x2l_ - dx_;
	//x2r_ = fsx2_ * scd_ / sad_;
	x2r_ = fsx2_;
	x1r_ = x2r_ + dx_;
	nl_.set(scd_, 0, x2l_);
	nr_.set(-scd_, 0, -x2r_);
	nl_.normalize();
	nr_.normalize();
	pl_.set(x2l_, 0, 0);
	pr_.set(x2r_, 0, 0);
}


double mcTransportJawPairFocused::getDistanceInside(mcParticle& p) const
{
	/*
	*/
	// HACK!! Из-за залипания частиц в сфокусированных объектах 
	// временно подменяем геометрию параллелепипедами
	geomVector3D pp(p.p);
	double x = pp.x();
	bool isRight = (x < fscenter_) ? false : true;

	if (isRight)
	{
		pp(0) = x - (fsx2_ + dx_ / 2);
		return mcGeometry::getDistanceToPrismInside(pp, p.u, dx_, dy_, h_);
	}
	else
	{
		pp(0) = x - (fsx1_ - dx_ / 2);
		return mcGeometry::getDistanceToPrismInside(pp, p.u, dx_, dy_, h_);
	}

	/*
	geomVector3D pp(p.p);
	geomVector3D v(p.u);

	double x = pp.x() * (scd_ / (scd_ - pp.z()));
	bool isRight = (x < fscenter_) ? false : true;

	double dx1 = DBL_MAX, dx2 = DBL_MAX;
	if (isRight)
	{
		double csn = v * nr_;
		if (csn > 0)
			dx2 = ((pr_ - pp) * nr_) / csn;
		if(v.x() > 0)
			dx1 = (x1r_ - pp.x()) / v.x();
	}
	else
	{
		double csn = v * nl_;
		if (csn > 0)
			dx2 = ((pl_ - pp) * nl_) / csn;
		if(v.x() < 0) 
			dx1 = (x1l_ - pp.x()) / v.x();
	}

	double dy = (v.y() < 0) ? fabs((pp.y() + 0.5 * dy_) / v.y()) : (v.y() > 0) ? fabs(((pp.y() - 0.5 * dy_)) / v.y()) : DBL_MAX;
	double dz = (v.z() < 0) ? fabs(pp.z() / v.z()) : (v.z() > 0) ? fabs((h_ - pp.z()) / v.z()) : DBL_MAX;
	return MIN(MIN(dx1, dx2), MIN(dy, dz));
	*/
}

double mcTransportJawPairFocused::getDistanceOutside(mcParticle& p) const
{
	/*
	*/
	// HACK!! Из-за залипания частиц в сфокусированных объектах 
	// временно подменяем геометрию параллелепипедами

	geomVector3D pp(p.p);
	double x = pp.x();

	pp(0) = x - (fsx1_ - dx_ / 2);
	double d1 = mcGeometry::getDistanceToPrismOutside(pp, p.u, dx_, dy_, h_);

	pp(0) = x - (fsx2_ + dx_ / 2);
	double d2 = mcGeometry::getDistanceToPrismOutside(pp, p.u, dx_, dy_, h_);

	return MIN(d1, d2);

	/*
	geomVector3D pp(p.p);
	geomVector3D v(p.u);
	double dist = 0;

	// Проверяем, не неаходится ли частица внутри описывающего параллелепипеда.
	double x = pp.x(), y = pp.y(), z = pp.z();

	bool isInside = x > x1l_ && x < x1r_ && fabs(y) < dy_ / 2 && z > 0 && z < h_;

	if (!isInside)
	{
		// Сначала дотаскиваем частицу до описывающего параллелепипеда
		double xc = (x1l_ + x1r_) / 2;
		geomVector3D pc(x - xc, y, z);
		double d = mcGeometry::getDistanceToPrismOutside(pc, v, x1r_ - x1l_, dy_, h_);
		if (d == DBL_MAX)
			return DBL_MAX;
		dist += d;
		pp += v * d;
		x = pp.x(); y = pp.y(); z = pp.z();

		// Проверяем, не попали ли в дырку
		if (v.z() > 0 && fabs(z) <= DBL_EPSILON) // частица летит вверх -> 
		{
			if (x > x2l_ && x < x2r_)
				isInside = true;
		}
		else if (v.z() < 0 && fabs(z - h_) <= DBL_EPSILON) // частица летит вниз -> 
		{
			if (x > (x2l_ * (scd_ - h_) / scd_) && x < (x2r_ * (scd_ - h_) / scd_))
				isInside = true;
		}
		else if ((v.y() < 0 && fabs(2 * y - dy_) <= DBL_EPSILON) ||
			(v.y() > 0 && fabs(2 * y + dy_) <= DBL_EPSILON))
		{
			double xx = x * (scd_ - x) / scd_;
			if (xx > x2l_&& xx < x2r_)
				isInside = true;
		}
	}

	if (isInside)
	{
		// Пересечение возможно только с внутренними гранями
		double csn = v * nl_;
		if (csn < 0)
		{
			double d = ((pl_ - pp) * nl_) / csn;
			pp += v * d;
			x = pp.x(); y = pp.y(); z = pp.z();
			if (fabs(pp.y()) < dy_ / 2 && pp.z() > 0 && pp.z() < h_)
				return dist + d;
			else
				return	DBL_MAX;
		}

		csn = v * nr_;
		if (csn < 0)
		{
			double d = ((pr_ - pp) * nr_) / csn;
			pp += v * d;
			x = pp.x(); y = pp.y(); z = pp.z();
			if (fabs(pp.y()) < dy_ / 2 && pp.z() > 0 && pp.z() < h_)
				return dist + d;
			else
				return	DBL_MAX;
		}

		// Летим от обеих граней
		return	DBL_MAX;
	}

	// Столкнулись с невнутренними гранями
	return dist;
	*/
}

double mcTransportJawPairFocused::getDNearInside(const geomVector3D& p) const
{
	return 0;
}

void mcTransportJawPairFocused::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "SAD:\t" << sad_ << endl;
	os << "SCD:\t" << scd_ << endl;
	os << "Dimensions:\t" << dx_ << "\t" << dy_ << endl;
	os << "Thickness:\t" << h_ << endl;
	os << "FSX:\t" << fsx1_ << "\t" << fsx2_ << endl;
}

void mcTransportJawPairFocused::dumpVRML(ostream& os) const
{
	os << "# Focused Jaws pair: " << this->getName() << endl;
	geomVector3D p[8];

	// HACK!! Симулируем рисование призм через отсутствие масштабирование верхней плоскости
	double f = 1;
	//double f = (scd_ - h_) / scd_;

	for (int k = 0; k < 2; k++)
	{
		int i = 0;
		if (k == 0)
		{
			// Левая часть
			p[i++] = geomVector3D(x1l_, -dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x1l_, dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x2l_, dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x2l_, -dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x1l_, -dy_ / 2, h_) * mttow_;
			p[i++] = geomVector3D(x1l_, dy_ / 2, h_) * mttow_;
			p[i++] = geomVector3D(x2l_ * f, dy_ / 2, h_) * mttow_;
			p[i++] = geomVector3D(x2l_ * f, -dy_ / 2, h_) * mttow_;
		}
		else
		{
			// Правая часть
			p[i++] = geomVector3D(x2r_, -dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x2r_, dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x1r_, dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x1r_, -dy_ / 2, 0) * mttow_;
			p[i++] = geomVector3D(x2r_ * f, -dy_ / 2, h_) * mttow_;
			p[i++] = geomVector3D(x2r_ * f, dy_ / 2, h_) * mttow_;
			p[i++] = geomVector3D(x1r_, dy_ / 2, h_) * mttow_;
			p[i++] = geomVector3D(x1r_, -dy_ / 2, h_) * mttow_;
		}

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

		for (i = 0; i < 8; i++) {
			os << "                    " << p[i].x() << ' ' << p[i].y() << ' ' << p[i].z();
			if (i < 7) os << ", ";
			os << endl;
		}

		os << "                ]" << endl;
		os << "            }" << endl;
		os << "            coordIndex [" << endl;

		os << "                0, 1, 2, 3, -1," << endl;
		os << "                0, 4, 5, 1, -1," << endl;
		os << "                1, 5, 6, 2, -1," << endl;
		os << "                2, 6, 7, 3, -1," << endl;
		os << "                3, 7, 4, 0, -1," << endl;
		os << "                4, 7, 6, 5, -1" << endl;

		os << "            ]" << endl;
		os << "        }" << endl;
		os << "      }" << endl;
		os << "    }" << endl;
	}
}
