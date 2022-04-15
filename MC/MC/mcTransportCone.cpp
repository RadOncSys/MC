#include "mcTransportCone.h"
#include <float.h>

mcTransportCone::mcTransportCone()
	:mcTransport()
{
	setGeometry(0, 0);
}

mcTransportCone::mcTransportCone(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double r, double h)
	: mcTransport(orgn, z, x)
{
	setGeometry(r, h);
}

mcTransportCone::~mcTransportCone(void)
{
}

void mcTransportCone::setGeometry(double r, double h)
{
	r_ = r;
	h_ = h;
	A_ = h == 0 ? 0 : r / h;
	A_ *= A_;
	g_.set(h_, r_);
	s_ = g_.length();
	g_.normalize();
}

double mcTransportCone::getDistanceInside(mcParticle& p)	const
{
	// Пересечение с основанием
	if (p.u.z() < 0) {
		double d = fabs(p.p.z() / p.u.z());
		geomVector3D pp = p.p + (p.u*d);
		if (pp.sqLengthXY() <= r_*r_)
			return d;
	}

	// Пересечение с основанием отсутствует. Ищем пересечение с конусом.
	double a = p.u.sqLengthXY() - A_*p.u.z()*p.u.z();
	double b = p.u.x()*p.p.x() + p.u.y()*p.p.y() - A_*(p.p.z() - h_)*p.u.z();
	double c = p.p.sqLengthXY() - A_*(p.p.z() - h_)*(p.p.z() - h_);
	double det = b*b - a*c;

	// Решение может отсутствовать только за пределами конуса,
	// что практически возможно в окрестности его поверхности.
	// Поэтому, возвращаем нулевое расстояние.
	if (det < 0) // при нулевом детерминанте попадаем точтно в вершину конуса
		return 0;

	// В этом месте b не может быть нулевым, так как это означало бы 
	// нулевой детерминант, обработанный выше
	if (a == 0)
	{
		if (b == 0) return 0;
		double d = -0.5 * c / b;
		if (d <= 0) return 0;
		else return d;
	}
	else
	{
		double d1 = (-b - sqrt(det)) / a, d2 = (-b + sqrt(det)) / a;
		if (d1 > 0 && d1 < d2) return d1;
		else return d2 < 0 ? 0 : d2;
	}
}

double mcTransportCone::getDistanceOutside(mcParticle& p) const
{
	if (p.p.z() <= 0)
	{
		if (p.u.z() <= 0)
			return DBL_MAX;
		double d = fabs(p.p.z() / p.u.z());
		geomVector3D pp = p.p + (p.u*d);
		if (pp.sqLengthXY() <= r_*r_)
			return d;
	}

	// Теперь точно знаем, что через торец проникновение невозможно.
	double a = p.u.sqLengthXY() - A_*p.u.z()*p.u.z();
	double b = p.u.x()*p.p.x() + p.u.y()*p.p.y() - A_*(p.p.z() - h_)*p.u.z();
	double c = p.p.sqLengthXY() - A_*(p.p.z() - h_)*(p.p.z() - h_);
	double det = b*b - a*c;

	// Нет пересечения с конусом
	if (det < 0)
		return DBL_MAX;

	double d1 = DBL_MAX, d2 = d1;

	if (a == 0)
		d1 = -c / (2 * b);
	else {
		if (a < 0) {
			d1 = (-b + sqrt(det)) / a;
			d2 = (-b - sqrt(det)) / a;
		}
		else {
			d1 = (-b - sqrt(det)) / a;
			d2 = (-b + sqrt(det)) / a;
		}
		if (d2 < 0) return DBL_MAX; // все пересечения сзади
	}

	// Попытка избежать залипания на поверхности.
	// Если частица на поверхности и летит в сторону, то пересечений не будет
	if (fabs(c) < 1e-14 && a < 0)
		return DBL_MAX;

	if (d1 > 0) {
		geomVector3D pp = p.p + (p.u*d1);
		if (pp.z() > 0 && pp.z() < h_)
			return d1;
	}

	if (d2 > 0 && d2 != DBL_MAX) {
		geomVector3D pp = p.p + (p.u*d2);
		if (pp.z() > 0 && pp.z() < h_)
			return d2;
	}

	return DBL_MAX;
}

double mcTransportCone::getDNearInside(const geomVector3D& p) const
{
	double rm = r_ * (h_ - p.z()) / h_;
	double r = p.lengthXY();
	double dr = fabs((rm - r)*g_.x());
	return MIN(dr, fabs(p.z() - h_));
}

void mcTransportCone::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Radius:\t" << getRadius() << endl;
	os << "Height:\t" << getHeight() << endl;
}

void mcTransportCone::dumpVRML(ostream& os) const
{
	geomVector3D p = geomVector3D(0, 0, h_*0.5) * mttow_;
	geomVector3D v(0, 0, 1);
	v(3) = 0;
	v = v * mttow_;

	os << "# Cone: " << this->getName() << endl;
	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	//os << "  rotation 1 0 0 1.5708" << endl;
	os << "  rotation 1 0 0 " << (1.5708*v.z()) << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Cone { " << endl;
	os << "                      bottomRadius " << r_ << endl;
	os << "                      height " << h_ << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
