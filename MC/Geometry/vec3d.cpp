#include "vec3d.h"
#include "text.h"

bool geomVector3D::isInsideBBox(const geomVector3D& b1, const geomVector3D& b2)const
{
	return p_[0] >= b1.x() && p_[0] <= b2.x() &&
		p_[1] >= b1.y() && p_[1] <= b2.y() &&
		p_[2] >= b1.z() && p_[2] <= b2.z();
}

istream& operator >> (istream& is, geomVector3D& v)
{
	string line;
	getline(is, line, '\n');
	GetFloatArray(line, v.p_, 4);
	return is;
}

ostream& operator << (ostream& os, const geomVector3D& v)
{
	os << v.p_[0] << '\t' << v.p_[1] << '\t' << v.p_[2] << '\t' << v.p_[3] << endl;
	return os;
}
