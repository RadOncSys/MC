#include "mcTransportRectangleTrap.h"

mcTransportRectangleTrap::mcTransportRectangleTrap(
	const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double width, double length)
	:mcTransport(orgn, z, x)
	, swidth_(width / 2)
	, slength_(length / 2)
{
}

mcTransportRectangleTrap::~mcTransportRectangleTrap(void)
{
}

void mcTransportRectangleTrap::beginTransport(mcParticle& p)
{
	if (p.u.z() <= 0)
		return;

	// ƒоталкиваем частицу до плоскости
	geomVector3D pp = p.p * mwtot_;
	geomVector3D uu = (p.p + p.u) * mwtot_;
	uu = uu - pp;
	double f = -pp.z() / uu.z();
	double x = pp.x() + uu.x() * f;
	if (fabs(x) > swidth_)
		return;
	double y = pp.y() + uu.y() * f;
	if (fabs(y) > slength_)
		return;

	if (nextTransport_ != nullptr)
		nextTransport_->beginTransport(p);
}

void mcTransportRectangleTrap::dumpVRML(ostream& os) const
{
	int i = 0;
	double r = 15;  // размера бокса
	geomVector3D p[8];

	p[i++] = geomVector3D(-r, -r, 0) * mttow_;
	p[i++] = geomVector3D(-r, r, 0) * mttow_;
	p[i++] = geomVector3D(r, r, 0) * mttow_;
	p[i++] = geomVector3D(r, -r, 0) * mttow_;
	p[i++] = geomVector3D(-swidth_, -slength_, 0) * mttow_;
	p[i++] = geomVector3D(-swidth_, slength_, 0) * mttow_;
	p[i++] = geomVector3D(swidth_, slength_, 0) * mttow_;
	p[i++] = geomVector3D(swidth_, -slength_, 0) * mttow_;

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

	os << "                0, 1, 5, 4, -1," << endl;
	os << "                1, 2, 6, 5, -1," << endl;
	os << "                2, 3, 7, 6, -1," << endl;
	os << "                3, 0, 4, 7, -1" << endl;

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}
