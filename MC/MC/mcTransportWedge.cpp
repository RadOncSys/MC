#include "mcTransportWedge.h"
#include <float.h>

mcTransportWedge::mcTransportWedge() : mcTransport()
{
	setGeometry(0, 0, 0);
}

mcTransportWedge::mcTransportWedge(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
	double ax, double ay, double az)
	: mcTransport(orgn, vz, vx)
{
	setGeometry(ax, ay, az);
}

mcTransportWedge::~mcTransportWedge(void)
{
}

void mcTransportWedge::setGeometry(double ax, double ay, double az)
{
	if (az <= 0)
		throw std::exception("Wedge height should be positive");
	ax_ = ax / 2;
	ay_ = ay / 2;
	az_ = az;
	n_.set(0, az, ay);
	n_.normalize();
	cy_.set(0, ay_, 0);
}

double mcTransportWedge::getDistanceInside(mcParticle& p) const
{
	// 5 плоскостей. Выбираем минимальное расстояние
	double dx0 = p.u.x() < 0 ? -(p.p.x() + ax_) / p.u.x() : DBL_MAX;
	double dx1 = p.u.x() > 0 ? (ax_ - p.p.x()) / p.u.x() : DBL_MAX;
	double dy0 = p.u.y() < 0 ? -(p.p.y() + ay_) / p.u.y() : DBL_MAX;
	double dz0 = p.u.z() < 0 ? -p.p.z() / p.u.z() : DBL_MAX;
	double cn = p.u * n_;
	double dz1 = cn > 0 ? ((cy_ - p.p) * n_) / cn : DBL_MAX;
	return MIN(MIN(MIN(dy0, dz0), MIN(dx0, dx1)), dz1);
}

double mcTransportWedge::getDistanceOutside(mcParticle& p) const
{
	// Перебираем все плоскости и определяем если есть столкновение то попадаем ли в рамку.
	// Если попадаем в рамку, то это уже сразу правильный вход.

	// наклонная грань
	double cn = p.u * n_;
	if (cn < 0)
	{
		double h = (p.p - cy_) * n_;
		if (h > 0)
		{
			double dist = -h / cn;
			geomVector3D pc = p.p + (p.u * dist);
			if (pc.x() > -ax_ && pc.x() < ax_ && pc.y() > -ay_ && pc.y() < ay_)
				return dist;
		}
	}

	// Z = const
	if (p.u.z() < 0 && p.p.z() > 0)
	{
		double dist = -p.p.z() / p.u.z();
		geomVector3D pc = p.p + (p.u * dist);
		if (pc.x() > -ax_ && pc.x() < ax_ && pc.y() > -ay_ && pc.y() < ay_)
			return dist;
	}

	// Y = const
	if (p.u.y() > 0 && p.p.y() <= -ay_)
	{
		double dist = -(p.p.y() + ay_) / p.u.y();
		geomVector3D pc = p.p + (p.u * dist);
		if (pc.x() > -ax_ && pc.x() < ax_ && pc.z() > 0 && pc.z() < az_)
			return dist;
	}

	// -X = const
	if (p.u.x() > 0 && p.p.x() <= -ax_)
	{
		double dist = -(p.p.x() + ax_) / p.u.x();
		geomVector3D pc = p.p + (p.u * dist);
		if (pc.y() > -ay_ && pc.y() < ay_ && pc.z() < az_ * (ay_ - p.p.y()) / (2 * ay_))
			return dist;
	}

	// +X = const
	if (p.u.x() < 0 && p.p.x() >= ax_)
	{
		double dist = (ax_ - p.p.x()) / p.u.x();
		geomVector3D pc = p.p + (p.u * dist);
		if (pc.y() > -ay_ && pc.y() < ay_ && pc.z() < az_ * (ay_ - p.p.y()) / (2 * ay_))
			return dist;
	}

	// Не попали ни в одну из граней
	return DBL_MAX;
}

double mcTransportWedge::getDNearInside(const geomVector3D& p) const
{
	return 0;
}

void mcTransportWedge::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "ax,ay,az:\t" << ax_ << '\t' << ay_ << '\t' << az_ << endl;
}

void mcTransportWedge::dumpVRML(ostream& os) const
{
	int i = 0;
	geomVector3D p[6];

	p[i++] = geomVector3D(-ax_, -ay_, 0) * mttow_;
	p[i++] = geomVector3D(-ax_, ay_, 0) * mttow_;
	p[i++] = geomVector3D(ax_, ay_, 0) * mttow_;
	p[i++] = geomVector3D(ax_, -ay_, 0) * mttow_;
	p[i++] = geomVector3D(-ax_, -ay_, az_) * mttow_;
	p[i++] = geomVector3D(ax_, -ay_, az_) * mttow_;

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

	for (i = 0; i < 6; i++) {
		os << "                    " << p[i].x() << ' ' << p[i].y() << ' ' << p[i].z();
		if (i < 5) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	os << "                0, 1, 2, 3, -1," << endl;
	os << "                0, 4, 1, -1," << endl;
	os << "                3, 2, 5, -1," << endl;
	os << "                0, 3, 5, 4, -1," << endl;
	os << "                1, 4, 5, 2, -1," << endl;

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;

}
