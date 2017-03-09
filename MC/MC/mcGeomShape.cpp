#include ".\mcgeomshape.h"
#include <float.h>

mcGeomShape::mcGeomShape()
	:nSides_(0)
	, Sides_(nullptr)
{
}

mcGeomShape::~mcGeomShape()
{
	for (unsigned int i = 0; i < nSides_; i++)
		delete Sides_[i];
	nSides_ = 0;
}

inline double mcGeomShape::getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const
{
	double d = DBL_MAX;
	for (unsigned int i = 0; i < nSides_; i++)
	{
		double di = Sides_[i]->getDistance(p, v, inside);
		if (di < 0) di = 0;
		if (d > di) d = di;
	}
	return d;
}

inline double mcGeomShape::getDistanceInside(const geomVector3D& p, const geomVector3D& v) const
{
	return getDistance(p, v, true);
}

inline double mcGeomShape::getDistanceOutside(const geomVector3D& p, const geomVector3D& v) const
{
	return getDistance(p, v, false);
}

inline double mcGeomShape::getDNear(const geomVector3D& p) const
{
	double d = DBL_MAX;
	for (unsigned int i = 0; i < nSides_; i++)
	{
		double di = Sides_[i]->getDNear(p);
		if (di < 0) di = 0;
		if (d > di) d = di;
	}
	return d;
}

inline double mcGeomShape::getDNearInside(const geomVector3D& p) const
{
	return getDNear(p);
}
