#include "mcETransportConvexPolygonCircle.h"
#include "mcGeometry.h"

mcETransportConvexPolygonCircle::mcETransportConvexPolygonCircle(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, std::vector<double>& pz, std::vector<double>& pr)
	: mcTransport(orgn, z, x)
	, pz_(pz)
	, pr_(pr)
{
}

double mcETransportConvexPolygonCircle::getDistanceInside(mcParticle& p) const
{
	double dist = mcGeometry::getDistanceToConvexPolygonCircleInside(p.p, p.u, pz_, pr_);

	if (internalTransport_ != nullptr)
	{
		// К сожалению частицу нужно переводить в систему координат объекта
		mcParticle pp(p);
		pp.p = pp.p * mtoe_;
		pp.u = p.u.transformDirection(mtoe_);
		double dist2 = internalTransport_->getDistanceOutside(pp);
		if (dist2 < dist)
		{
			p.transportNearest_ = pp.transportNearest_;
			p.exitSurface_ = mcParticle::temb_shit_t::Internal;
			return dist2;
		}
		else
		{
			p.exitSurface_ = mcParticle::temb_shit_t::External;
			return dist;
		}
	}
	else
	{
		p.exitSurface_ = mcParticle::temb_shit_t::External;
		return dist;
	}
}

double mcETransportConvexPolygonCircle::getDistanceOutside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;
	return mcGeometry::getDistanceToConvexPolygonCircleOutside(p.p, p.u, pz_, pr_);
}

void mcETransportConvexPolygonCircle::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Polygon:\t" << endl;
	unsigned i;
	for (i = 0; i < pz_.size(); i++) os << pz_[i] << "\t";
	os << endl;
	for (i = 0; i < pr_.size(); i++) os << pr_[i] << "\t";
	os << endl << endl;
}

void mcETransportConvexPolygonCircle::dumpVRML(ostream& os) const
{
	os << "# mcETransportConvexPolygonCircle: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	dumpVRMLPolygonCircle(os, pz_, pr_);

	os << "  ]" << endl;
	os << "}" << endl;
}
