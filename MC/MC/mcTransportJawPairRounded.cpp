#include "mcTransportJawPairRounded.h"
#include "mcGeometry.h"
#include <float.h>

mcTransportJawPairRounded::mcTransportJawPairRounded(void)
	:mcTransport()
{
}

mcTransportJawPairRounded::mcTransportJawPairRounded(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x,
	double r, double h, double dx)
	: mcTransport(orgn, z, x)
	, r_(r), h_(h), dx_(dx), r2_(r*r)
{
}

mcTransportJawPairRounded::~mcTransportJawPairRounded(void)
{
}

double mcTransportJawPairRounded::getDistanceInside(mcParticle& p) const
{
	double x = p.p.x();
	bool isRight = (x < fscenter_) ? false : true;
	double dist = 0;

	// Приводим расчет только к одной половинке
	geomVector3D position(p.p);
	geomVector3D direction(p.u);
	if (isRight)
	{
		position(0) = fsx1_ + fsx2_ - position.x();
		position(1) = -position.y();
		direction(0) = -direction.x();
		direction(1) = -direction.y();
		x = position.x();
	}

	// Нужно двигаться в два этапа.
	// Сначала на поверхность одной части, затем, в зависимости от места пересечения, добавляем отрезок внутри второй.
	x -= fsx1_ - dx_;
	geomVector3D pp(position);
	pp(0) = x;
	if (x < 0)
	{
		double d = mcGeometry::getDistanceToCylinderInside(pp, direction, r_, h_);
		geomVector3D p1 = pp + (direction * d);
		if (p1.x() >= 0)
		{
			dist = d * x / (x - p1.x());
			p1 = pp + (direction * dist);
			p1(0) -= dx_ / 2;
			d = mcGeometry::getDistanceToPrismInside(p1, direction, dx_, r_ * 2, h_);
			dist += d;
		}
		else
			dist = d;
	}
	else
	{
		pp(0) = x - dx_ / 2;
		dist = mcGeometry::getDistanceToPrismInside(pp, direction, dx_, r_ * 2, h_);
		geomVector3D p1 = pp + (direction * dist);
		p1(0) += dx_ / 2;
		if (p1.x() <= 0)
		{
			dist += mcGeometry::getDistanceToCylinderInside(p1, direction, r_, h_);
		}
	}

	return dist;
}

double mcTransportJawPairRounded::getDistanceOutside(mcParticle& p) const
{
	double cd1 = getDistanceToLeftOutside(p.p, p.u);

	geomVector3D pp(fsx1_ + fsx2_ - p.p.x(), -p.p.y(), p.p.z());
	geomVector3D uu(-p.u.x(), -p.u.y(), p.u.z());
	double cd2 = getDistanceToLeftOutside(pp, uu);

	return MIN(cd1, cd2);
}

double mcTransportJawPairRounded::getDistanceToLeftOutside(const geomVector3D& p, const geomVector3D& u) const
{
	geomVector3D pp(p);
	pp(0) = p.x() - fsx1_ + dx_;
	double cd1 = mcGeometry::getDistanceToCylinderOutside(pp, u, r_, h_);
	pp(0) -= dx_ / 2;
	double cd2 = mcGeometry::getDistanceToPrismOutside(pp, u, dx_, r_ * 2, h_);
	if (cd1 == DBL_MAX || p.x() + cd1 * u.x() > fsx1_) 
		return cd2;
	else 
		return MIN(cd1, cd2);
}

double mcTransportJawPairRounded::getDNearInside(const geomVector3D& p) const
{
	return 0;
}

void mcTransportJawPairRounded::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Radius step:\t" << r_ << endl;
	os << "Thickness:\t" << h_ << endl;
	os << "Rectangle part width:\t" << dx_ << endl;
}

void mcTransportJawPairRounded::dumpVRML(ostream& os) const
{
	os << "# Cylindrical grid: " << this->getName() << endl;

	for (int k = 0; k < 2; k++)
	{
		//geomMatrix3D M = (k == 0) ? (mttow_ * geomMatrix3D::ParallelShift(fsx1_ - dx_, 0, 0)) :
		//	(mttow_ * geomMatrix3D::ParallelShift(fsx2_ + dx_, 0, 0) * geomMatrix3D::RotationAroundOz(180));
		geomMatrix3D M = (k == 0) ? (geomMatrix3D::ParallelShift(fsx1_ - dx_, 0, 0) * mttow_) :
			(geomMatrix3D::RotationAroundOz(180) * geomMatrix3D::ParallelShift(fsx2_ + dx_, 0, 0) * mttow_);

		// Прямоугольная часть
		int i = 0;
		geomVector3D p[8];
		p[i++] = geomVector3D(0, -r_, 0) * M;
		p[i++] = geomVector3D(0, r_, 0) * M;
		p[i++] = geomVector3D(dx_, r_, 0) * M;
		p[i++] = geomVector3D(dx_, -r_, 0) * M;
		p[i++] = geomVector3D(0, -r_, h_) * M;
		p[i++] = geomVector3D(0, r_, h_) * M;
		p[i++] = geomVector3D(dx_, r_, h_) * M;
		p[i++] = geomVector3D(dx_, -r_, h_) * M;

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
		os << "                1, 5, 6, 2, -1," << endl;
		os << "                2, 6, 7, 3, -1," << endl;
		os << "                3, 7, 4, 0, -1," << endl;
		os << "                4, 7, 6, 5, -1" << endl;

		os << "            ]" << endl;
		os << "        }" << endl;
		os << "      }" << endl;
		os << "    }" << endl;

		// Полуцилиндр
		dumpVRMLCylinderSemiSide(os, r_, h_, M);
		dumpVRMLSemiCircle(os, r_, h_, true, M);
		dumpVRMLSemiCircle(os, r_, 0, false, M);
	}
}
