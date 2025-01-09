#include "mcgeomrectside.h"
#include "mcDefs.h"
#include <float.h>

mcGeomRectSide::mcGeomRectSide()
{
}

mcGeomRectSide::mcGeomRectSide(const geomVector3D& p, 
	const geomVector3D& Vx, const geomVector3D& Vy, double ax, double ay)
{
	SetGeometry(p, Vx, Vy, ax, ay);
}

void mcGeomRectSide::SetGeometry(const geomVector3D& p,
	const geomVector3D& Vx, const geomVector3D& Vy, double ax, double ay)
{
	ax_ = ax / 2; ay_ = ay / 2;

	geomVector3D vx(Vx);
	geomVector3D vy(Vy);
	geomVector3D vz = Vx ^ Vy;
	vx.normalize();
	vy.normalize();
	vz.normalize();

	m_ = geomMatrix3D::ParallelShift(-p.x(), -p.y(), -p.z()) *
		geomMatrix3D::BuildFromAxis(vx, vy, vz);
}

double mcGeomRectSide::getDistance(const geomVector3D& p, const geomVector3D& v) const
{
	auto pp = p * m_;
	auto vv = v.transformDirection(m_);

	// От плоскости
	if ((pp.z() * vv.z()) >= 0)
		return DBL_MAX;
	else
	{
		double dist = -pp.z() / vv.z();
		auto s = pp + (vv * dist);

		// Попадание в рамку
		if (ABS(s.x()) <= ax_ && ABS(s.y()) <= ay_)
			return dist;
		else
			return DBL_MAX;
	}
}

double mcGeomRectSide::getDNear(const geomVector3D& p) const
{
	auto pp = p * m_;
		
	if (ABS(pp.x()) <= ax_ && ABS(pp.y()) <= ay_)
		return ABS(pp.z());
	else if (ABS(pp.x()) <= ax_)
	{
		float dy = ABS(pp.y()) - ay_;
		return sqrt(pp.z() * pp.z() + dy * dy);
	}
	else if (ABS(pp.y()) <= ay_)
	{
		float dx = ABS(pp.x()) - ax_;
		return sqrt(pp.z() * pp.z() + dx * dx);
	}
	else
	{
		float dx = ABS(pp.x()) - ax_;
		float dy = ABS(pp.y()) - ay_;
		return sqrt(pp.z() * pp.z() + dx * dx + dy * dy);
	}
}
