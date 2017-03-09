#include "mcTransportPrism.h"
#include "mcGeometry.h"

mcTransportPrism::mcTransportPrism()
	:mcTransport()
{
	setGeometry(0, 0, 0);
}

mcTransportPrism::mcTransportPrism(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
	double ax, double ay, double az)
	: mcTransport(orgn, vz, vx)
{
	setGeometry(ax, ay, az);
}

mcTransportPrism::~mcTransportPrism(void)
{
}

void mcTransportPrism::setGeometry(double ax, double ay, double az)
{
	ax_ = ax;
	ay_ = ay;
	az_ = az;
}

double mcTransportPrism::getDistanceInside(mcParticle& p) const
{
	return mcGeometry::getDistanceToPrismInside(p.p, p.u, ax_, ay_, az_);
}

double mcTransportPrism::getDistanceOutside(mcParticle& p) const
{
	return mcGeometry::getDistanceToPrismOutside(p.p, p.u, ax_, ay_, az_);
}

double mcTransportPrism::getDNearInside(const geomVector3D& p) const
{
	// HACK!!
	// ѕо непон€тным причинам использование DNEAR приводит к проблемам вблизи поверхностей.
	// Ёто серьезный предмет дл€ разбирательства, так как код ниже, похоже работает правильно
	// и причина где-то глубоко внутри движка транспорта.
	// ѕоследнее чревато фатальными проблемами.
	return 0;

	double dnear = MIN(fabs(p.z()), fabs(p.z() - az_));
	dnear = MIN(fabs(p.x() - 0.5*ax_), dnear);
	dnear = MIN(fabs(p.x() + 0.5*ax_), dnear);
	dnear = MIN(fabs(p.y() - 0.5*ay_), dnear);
	dnear = MIN(fabs(p.y() + 0.5*ay_), dnear);
	return dnear;
}

void mcTransportPrism::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "ax,ay,az:\t" << ax_ << '\t' << ay_ << '\t' << az_ << endl;
}

void mcTransportPrism::dumpVRML(ostream& os) const
{
	dumpVRMLPrism(os, ax_, ay_, az_);
}
