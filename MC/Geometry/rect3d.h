// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "frect.h"
#include "plane3d.h"

class geomRect3D :public geomPlane3D, public geomFRect
{
public:
	geomRect3D();

	geomRect3D(const geomVector3D& p
		, const geomVector3D& n
		, const geomVector3D& xv
		, const geomFRect& rect) {
		set(p, n, xv, rect);
	}

	geomRect3D(const geomRect3D& r);

	void set(const geomVector3D&
		, const geomVector3D&
		, const geomVector3D&
		, const geomFRect&);

	geomVector3D getCenter()const;

	void operator=(const geomRect3D&);

	friend istream& operator >> (istream&, geomRect3D&);
	friend ostream& operator << (ostream&, const geomRect3D&);
};
