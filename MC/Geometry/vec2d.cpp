#include "vec2d.h"

geomVector2D
crossTwo2DLines(const geomVector2D& p1
	, const geomVector2D& p2
	, const geomVector2D& p3
	, const geomVector2D& p4)
{
	geomVector2D n = p2 - p1;
	n.turnRight();
	double h0 = (p3 - p1)*n, h1 = (p4 - p1)*n;
	if (h0 == h1)
		throw std::exception("crossTwo2DLines: линии параллельны");
	return p3 + (p4 - p3)*(h0 / (h0 - h1));
}
