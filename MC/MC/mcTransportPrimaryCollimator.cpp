#include "mcTransportPrimaryCollimator.h"
#include <float.h>

mcTransportPrimaryCollimator::mcTransportPrimaryCollimator(void)
{
	setGeometry(0, 0, 0);
}
mcTransportPrimaryCollimator::mcTransportPrimaryCollimator(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x,
	double r0, double r1, double h) :mcTransport(orgn, z, x)
{
	setGeometry(r0, r1, h);
}

mcTransportPrimaryCollimator::~mcTransportPrimaryCollimator(void) {}

void mcTransportPrimaryCollimator::setGeometry(double r0, double r1, double h)
{
	r0_ = r0;
	r1_ = r1;
	h_col_ = h;//толщина коллиматора
	A_ = h == 0 ? 0 : (r1 - r0) / h;
	A_ *= A_;
}

double mcTransportPrimaryCollimator::getDistanceInside(mcParticle& p) const
{
	/*int static count01=0;
	count01++;
	if(count01==1)
	{
	cout << "InsideMethod" << endl << endl;
	}*/

	double px = p.p.x(), py = p.p.y(), pz = p.p.z();
	double vx = p.u.x(), vy = p.u.y(), vz = p.u.z();

	double dCone, dSlab;
	static double h = h_col_*(r1_ / (r1_ - r0_));//высота конуса вращения

	double a = vx*vx + vy*vy - A_*vz*vz;
	double b = px*vx + py*vy - A_*vz*(pz + (h - h_col_));
	double c = px*px + py*py - A_*(pz + (h - h_col_))*(pz + (h - h_col_));
	double det = b*b - a*c;

	dSlab = (vz < 0) ? (-pz) / vz : (h_col_ - pz) / vz;

	if ((a == 0) || (det <= 0))//ч-ца внутри кол-ра летит по касательной к конусу
		dCone = DBL_MAX;
	else
		dCone = (a > 0) ? ((-b - sqrt(det)) / a) : ((-b + sqrt(det)) / a);

	dCone = (dCone < 0) ? DBL_MAX : dCone;

	return MIN(dCone, dSlab);
}

double mcTransportPrimaryCollimator::getDistanceOutside(mcParticle& p) const
{
	/*int static count01=0;
	count01++;
	if(count01==1)
	{
	cout << "InsideMethod" << endl << endl;
	}*/

	double px = p.p.x(), py = p.p.y(), pz = p.p.z();
	double vx = p.u.x(), vy = p.u.y(), vz = p.u.z();

	if ((pz <= 0 && vz <= 0) || (pz >= h_col_ && vz >= 0))
		return DBL_MAX;


	double dCone1, dCone2, dSlab;
	static double h = h_col_*(r1_ / (r1_ - r0_));//высота конуса вращения

	double a = vx*vx + vy*vy - A_*vz*vz;
	double b = px*vx + py*vy - A_*vz*(pz + (h - h_col_));
	double c = px*px + py*py - A_*(pz + (h - h_col_))*(pz + (h - h_col_));
	double det = b*b - a*c;

	if (pz<0 || pz>h_col_)
		dSlab = (pz < 0) ? (-pz) / vz : (h_col_ - pz) / vz;
	else
		dSlab = DBL_MAX;

	if ((a == 0) || (det <= 0))
		dCone1 = dCone2 = DBL_MAX;
	else
	{
		dCone1 = (a > 0) ? (-b - sqrt(det)) / a : (-b + sqrt(det)) / a;
		dCone2 = (a > 0) ? (-b + sqrt(det)) / a : (-b - sqrt(det)) / a;
	}

	geomVector3D pp_Slab = p.p + (p.u*dSlab);
	geomVector3D pp_Cone1 = p.p + (p.u*dCone1);
	geomVector3D pp_Cone2 = p.p + (p.u*dCone2);

	if (((pp_Slab.z() < 0.1) && (pp_Slab.sqLengthXY() > r0_*r0_)) ||//pp_Slab.z() < 0.1 ----> избежание залипания частицы на поверхности
		((pp_Slab.z() > (h_col_ - 0.1)) && (pp_Slab.sqLengthXY() > r1_*r1_)))//pp_Slab.z() > (h_col_-0.1) ----> избежание залипания частицы на поверхности
		dSlab = dSlab;
	else
		dSlab = DBL_MAX;


	dCone1 = ((pp_Cone1.z() < 0) || (pp_Cone1.z() > h_col_) || (dCone1 < 0)) ? DBL_MAX : dCone1;
	dCone2 = ((pp_Cone2.z() < 0) || (pp_Cone2.z() > h_col_) || (dCone2 < 0)) ? DBL_MAX : dCone2;

	//избежание залипания частицы на поверхности
	if (fabs(c) < 1e-10 && a*b < 0)
		dCone1 = DBL_MAX;

	dCone1 = MIN(dCone1, dCone2);

	return MIN(dCone1, dSlab);
}

//return later
double mcTransportPrimaryCollimator::getDNearInside(const geomVector3D& p) const
{
	return 0;
}

void mcTransportPrimaryCollimator::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Minor Radius:\t" << R0() << endl;
	os << "Major Radius:\t" << R1() << endl;
	os << "Height(Thickness):\t" << getHeight() << endl;
}

void mcTransportPrimaryCollimator::dumpVRML(ostream& os) const
{
	double h_ = h_col_ * (r1_ / (r1_ - r0_));//высота конуса
	geomVector3D p = geomVector3D(0, 0, h_*0.5) * mttow_;
	geomVector3D v(0, 0, 1);
	v(3) = 0;
	v = v * mttow_;

	os << "# Collimator: " << this->getName() << endl;
	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	os << "  rotation 1 0 0 " << (1.5708*v.z()) << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Collimator { " << endl;
	os << "                      Minor Radius " << r0_ << endl;
	os << "                      Major Radius " << r1_ << endl;
	os << "                      height " << h_ << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}