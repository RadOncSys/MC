#include "mcTransportPolygon.h"
#include "../geometry/vec2d.h"

using namespace std;

mcTransportPolygon::mcTransportPolygon(void) : _ptr(nullptr), _np(0), _r(0), _h(0)
{
}


mcTransportPolygon::~mcTransportPolygon(void)
{
}

mcTransportPolygon::mcTransportPolygon(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, int np, double h, double r)
{
	setGeometry(np, h, r);
}

//return later
void mcTransportPolygon::setGeometry(int np, double h, double r)
{
	_np = np;	_h = h;		_r = r;

	geomVector2D* polyArray = new geomVector2D[_np];
	_ptr = polyArray;

	polyArray[0] = geomVector2D(h, 0);
	polyArray[1] = geomVector2D(4, 2);
	polyArray[2] = geomVector2D(2, 4);
	polyArray[3] = geomVector2D(0, r);

	//cout << "========================================================";
	//cout << "Insert (z,r)-coordinates of polygon : " << endl;

	//for(int i = 0; i<_np; i++)
	//{
	//	if(i==0)
	//	{
	//		cout << "top point ";
	//		double zTop, rTop;
	//		cin >> zTop >> rTop;
	//		polyArray[i] = geomVector2D(zTop, rTop);//?проверить, что пишется в массив
	//	}
	//	else
	//	{
	//		cout << (i+1) << " point ";
	//		double z_, r_;
	//		cin >> z_ >> r_;
	//		polyArray[i] = geomVector2D(z_, r_);
	//	}
	//}

	//cout << "========================================================";

	//проверить как выполняется эта строка??
	//_r = _ptr->y();
}

double mcTransportPolygon::getDistanceInside(mcParticle& p) const
{
	double vx = p.u.x(), vy = p.u.y(), vz = p.u.z();
	double px = p.p.x(), py = p.p.y(), pz = p.p.z();

	if (vz < 0)//пересечение с основанием 
	{
		double d = fabs(pz / vz);
		geomVector3D pp = p.p + (p.u*d);
		if (pp.sqLengthXY() <= _r*_r)
			return d;
	}

	double dCone = DBL_MAX;
	double dCone_temp;

	//пересечение с основанием отсутствует;
	//ищем пересечение с боковой поверхностью.
	for (int i = 0; i < (_np - 1); i++)
	{
		geomVector2D* ptr = _ptr;
		double z_i = ptr->x();
		double r_i = ptr->y();
		ptr++;
		double z_i_plus1 = ptr->x();
		double r_i_plus1 = ptr->y();
		double A = (r_i - r_i_plus1) / (z_i - z_i_plus1);//тангенс угла наклона одного из конусов, образующих полигон 
		A *= A;

		double a = p.u.sqLengthXY() - A*vz*vz;
		double b = px*vx + py*vy - A*vz*(pz - _h);
		double c = p.p.sqLengthXY() - A*(pz - _h)*(pz - _h);
		double det = b*b - a*c;

		// Решение может отсутствовать только за пределами конуса,
		// что практически возможно в окрестности его поверхности.
		if (det < 0)
			dCone_temp = DBL_MAX;
		//пересечение через вершину
		else if (det == 0)
			dCone_temp = -b / a;
		//пересечение через боковую поверхность
		else
			dCone_temp = (a > 0) ? ((-b + sqrt(det)) / a) : ((-b - sqrt(det)) / a);

		dCone = (dCone > dCone_temp) ? dCone_temp : dCone;
	}

	return dCone;
}


double mcTransportPolygon::getDistanceOutside(mcParticle& p) const
{
	double vx = p.u.x(), vy = p.u.y(), vz = p.u.z();
	double px = p.p.x(), py = p.p.y(), pz = p.p.z();

	//частица летит от полигона
	if ((pz <= 0 && vz <= 0) || (pz >= _h && vz >= 0))
		return DBL_MAX;

	//пересечение основаниe полигона
	if (pz <= 0)
	{
		double d = fabs(pz / vz);
		geomVector3D pp = p.p + (p.u*d);
		if (pp.sqLengthXY() <= _r*_r)
			return d;
	}

	//пересечение с боковой поверхностью полигона
	double dCone = 0;
	double dCone_temp;

	for (int i = 0; i < (_np - 1); i++)
	{
		geomVector2D* ptr = _ptr;
		double z_i = ptr->x();
		double r_i = ptr->y();
		ptr++;
		double z_i_plus1 = ptr->x();
		double r_i_plus1 = ptr->y();
		double A = (r_i - r_i_plus1) / (z_i - z_i_plus1);//тангенс угла наклона одного из конусов, образующих полигон 
		A *= A;

		double a = p.u.sqLengthXY() - A*vz*vz;
		double b = px*vx + py*vy - A*vz*(pz - _h);
		double c = p.p.sqLengthXY() - A*(pz - _h)*(pz - _h);
		double det = b*b - a*c;

		//пересечений нет
		if (det <= 0 || a == 0)
			dCone = 0;
		//пересечение через боковую поверхность
		else
		{
			dCone_temp = (a > 0) ? (-b - sqrt(det)) / a : (-b + sqrt(det)) / a;
			geomVector3D pp = p.p + (p.u*dCone_temp);
			if ((pp.z() >= z_i_plus1) && (pp.z() <= z_i))
				dCone = dCone_temp;
		}
	}

	return dCone;
}

//return later
double mcTransportPolygon::getDNearInside(const geomVector3D& p) const
{
	return 0;
}
