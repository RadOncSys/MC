#include "mcGeomWedgeTriangle.h"
#include "mcGeomRectSide.h"
#include "mcGeomTriangleSide.h"

mcGeomWedgeTriangle::
mcGeomWedgeTriangle(double x0, double y0, double x1, double xw)
	:x0_(x0)
	, y0_(y0)
	, x1_(x1)
	, xw_(xw / 2)
{
	//setObjName("Wedge");

	// Грани
	nSides_ = 5;
	Sides_ = new mcGeomSide*[nSides_];

	Sides_[0] = new mcGeomRectSide(geomVector3D(-xw_, x0_, 0)
		, geomVector3D(0, 1, 0)
		, geomVector3D(1, 0, 0)
		, x1_ - x0_, xw_ * 2);

	Sides_[1] = new mcGeomRectSide(geomVector3D(-xw_, x0_, 0)
		, geomVector3D(1, 0, 0)
		, geomVector3D(0, 0, 1)
		, xw_ * 2, y0_);

	geomVector3D tv = geomVector3D(0, x0_ - x1_, y0_);
	Sides_[2] = new mcGeomRectSide(geomVector3D(xw_, x1_, 0)
		, geomVector3D(-1, 0, 0)
		, tv
		, xw_ * 2, tv.length());

	Sides_[3] = new mcGeomTriangleSide(geomVector3D(-xw_, x0_, 0)
		, geomVector3D(0, 0, 1)
		, geomVector3D(0, 1, 0)
		, y0_, x1_ - x0_);

	Sides_[4] = new mcGeomTriangleSide(geomVector3D(xw_, x0_, 0)
		, geomVector3D(0, 1, 0)
		, geomVector3D(0, 0, 1)
		, x1_ - x0_, y0_);
}
