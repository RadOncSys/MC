#include "grid3d.h"
#include <string>

#ifndef ROUND
#define	ROUND(x)		(((x)>=0.0)?int((x)+0.5):int((x)-0.5))
#endif

geomGrid3D::geomGrid3D()
	:nx_(0)
	, ny_(0)
	, psx_(0)
	, psy_(0)
{}

geomGrid3D::
geomGrid3D(const geomGrid3D& g)
	:geomRect3D(g)
	, nx_(g.numCol())
	, ny_(g.numRow())
	, psx_(g.pixelSX())
	, psy_(g.pixelSY())
{}

geomVector3D
geomGrid3D::
getPointByIndex(int col, int row)const
{
	double x = minX() + (0.5 + double(col))*psx_;
	double y = minY() + (0.5 + double(row))*psy_;
	return geomVector3D(x, y, 0)*mRtoP_;
}

void
geomGrid3D::
getPointByIndex(int col, int row, double& x, double& y)const
{
	x = minX() + (0.5 + double(col))*psx_;
	y = minY() + (0.5 + double(row))*psy_;
}

void
geomGrid3D::
setGrid(int nx, int ny, double psx, double psy) {
	nx_ = nx;
	ny_ = ny;
	setResolution(psx, psy);
}

void
geomGrid3D::
setResolution(double x, double y) {
	psx_ = x;
	psy_ = y;
	if (r_ == l_) { r_ = 0.5*nx_*psx_; l_ = -r_; }
	if (t_ == b_) { t_ = 0.5*ny_*psy_; b_ = -t_; }
}

void
geomGrid3D::
setAdjustedResolution(double x, double y) {
	psx_ = x;
	psy_ = y;
	l_ = x*(ROUND(l_ / x) - 0.5);
	b_ = y*(ROUND(b_ / y) - 0.5);
	r_ = x*(ROUND(r_ / x) + 0.5);
	t_ = y*(ROUND(t_ / y) + 0.5);
	nx_ = ROUND((r_ - l_) / x);
	ny_ = ROUND((t_ - b_) / y);
}

void
geomGrid3D::
setSize(int nx, int ny) {
	nx_ = nx; ny_ = ny; psx_ = width() / nx_; psy_ = height() / ny_;
}

int
geomGrid3D::
getColIndex(double x)const
{
	return int(nx_*(x - l_) / (r_ - l_));
}

int
geomGrid3D::
getRowIndex(double y)const
{
	return int(ny_*(y - b_) / (t_ - b_));
}

double
geomGrid3D::
getCol(unsigned i)const
{
	return minX() + (0.5 + double(i))*psx_;
}

double
geomGrid3D::
getRow(unsigned j)const
{
	return minY() + (0.5 + double(j))*psy_;
}

void geomGrid3D::scale(double f)
{
	geomFRect::scale(f);
	psx_ *= f;
	psy_ *= f;
}

void geomGrid3D::operator=(const geomGrid3D& g)
{
	*(geomRect3D*)this = (const geomRect3D&)g;
	nx_ = g.nx_;
	ny_ = g.ny_;
	psx_ = g.psx_;
	psy_ = g.psy_;
}

istream& operator >> (istream& is, geomGrid3D& g)
{
	is >> (geomRect3D&)g;
	string line;
	getline(is, line, '\n'); g.nx_ = atoi(line.c_str());
	getline(is, line, '\n'); g.ny_ = atoi(line.c_str());
	getline(is, line, '\n'); g.psx_ = atof(line.c_str());
	getline(is, line, '\n'); g.psy_ = atof(line.c_str());
	return is;
}

ostream& operator << (ostream& os, const geomGrid3D& g)
{
	os << (geomRect3D&)g;
	os << g.nx_ << endl;
	os << g.ny_ << endl;
	os << g.psx_ << endl;
	os << g.psy_ << endl;
	return os;
}
