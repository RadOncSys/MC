#include "mcPTLasVegas.h"
#include "mcGeometry.h"
#include <float.h>

mcPTLasVegas::mcPTLasVegas(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x) : mcTransport(orgn, z, x)
, a_(14.0)
, z_(2.175)
, ds_(6)
, hs_(5)
, xs_(6)
, ys_(5)
, x1_(0), x2_(0), y1_(0), y2_(0), z2_(0)
{
	ds_[0] = 0.1; ds_[1] = 0.2; ds_[2] = 0.4; ds_[3] = 0.7; ds_[4] = 1.0; ds_[5] = 1.5;
	hs_[0] = 0.05; hs_[1] = 0.1; hs_[2] = 0.2; hs_[3] = 0.325; hs_[4] = 0.45;
	xs_[0] = -4.5; xs_[1] = -3.0; xs_[2] = -1.5; xs_[3] = 0; xs_[4] = 2.0; xs_[5] = 4.0;
	ys_[0] = -4.0; ys_[1] = -2.0; ys_[2] = 0; ys_[3] = 2.0; ys_[4] = 4.0;

	unsigned i, j;
	for (j = 0; j < ys_.size(); j++)
	{
		if (hs_[j] > z2_) z2_ = hs_[j];
		for (i = 0; i < xs_.size(); i++)
		{
			if (x1_ > xs_[i] - ds_[i]) x1_ = xs_[i] - ds_[i];
			if (x2_ < xs_[i] + ds_[i]) x2_ = xs_[i] + ds_[i];
			if (y1_ > ys_[j] - ds_[i]) y1_ = ys_[j] - ds_[i];
			if (y2_ < ys_[j] + ds_[i]) y2_ = ys_[j] + ds_[i];
		}
	}
	ax_ = 2.0 * MAX(-x1_, x2_);
	ay_ = 2.0 * MAX(-y1_, y2_);

	//// “ест со значительно большими толщинами
	//for (i = 0; i < hs_.size(); i++) hs_[i] *= 2.2;
}

mcPTLasVegas::~mcPTLasVegas(void)
{
}

double mcPTLasVegas::getDistanceInside(mcParticle& p) const
{
	// –ассто€ние до параллелепипеда нужно в любом случае
	double dprism = mcGeometry::getDistanceToPrismInside(p.p, p.u, a_, a_, z_);

	// ƒл€ ускорени€ определ€ем пересечение или нахождение в области отверстий.
	double x = p.p.x(), y = p.p.y(), z = p.p.z();
	if ((fabs(x) <= ax_ / 2 && fabs(y) <= ay_ / 2 && fabs(z) <= z2_) ||
		mcGeometry::getDistanceToPrismOutside(p.p, p.u, ax_, ay_, z2_) != DBL_MAX)
	{
		geomVector3D pp = p.p;

		// ѕровер€ем пересечение с цилиндрами
		double dcyl = dprism;
		unsigned i, j;
		for (j = 0; j < ys_.size(); j++)
		{
			pp(1) = p.p.y() - ys_[j];
			for (i = 0; i < xs_.size(); i++)
			{
				pp(0) = p.p.x() - xs_[i];
				double d = mcGeometry::getDistanceToCylinderOutside(pp, p.u, ds_[i] / 2, hs_[j]);
				if (d < 0)
					d = mcGeometry::getDistanceToCylinderOutside(pp, p.u, ds_[i] / 2, hs_[j]);
				if (dcyl > d) dcyl = d;
			}
		}
		return dcyl;
	}
	else
		return dprism;
}

double mcPTLasVegas::getDistanceOutside(mcParticle& p) const
{
	// —начала определе€ем пересечение с параллелепипедом фантом и перемещаем частицу на его поверхность.
	// «атем определ€емс€, не неаходитс€ ли она в рамке, где могут быть высверленные отверсти€.

	// Ќа случай, если частица внутри паралелепида в каком-то отверстии
	double x = p.p.x(), y = p.p.y(), z = p.p.z();
	bool ispinside = fabs(2 * x) < a_ - FLT_EPSILON && fabs(2 * y) < a_ - FLT_EPSILON && z > FLT_EPSILON && z < z_ / 2;

	if (ispinside)
	{
		unsigned i, j;
		// »щем цилиндр
		for (i = 0; i < xs_.size(); i++)
			if (fabs(x - xs_[i]) <= ds_[i] / 2) break;
		if (i >= xs_.size())
			throw std::exception("mcPTLasVegas::getDistanceOutside: Particle expected to be inside cylinder, but cylinder not found");

		for (j = 0; j < ys_.size(); j++)
			if (((x - xs_[i])*(x - xs_[i]) + (y - ys_[j])*(y - ys_[j]) - FLT_EPSILON) * 4 <= ds_[i] * ds_[i]) break;
		if (j >= ys_.size())
			throw std::exception("mcPTLasVegas::getDistanceOutside: Particle expected to be inside cylinder, but cylinder not found");

		geomVector3D pp = p.p;
		pp(0) -= xs_[i]; pp(1) -= ys_[j];
		double dcyl = mcGeometry::getDistanceToCylinderInside(pp, p.u, ds_[i] / 2, hs_[j]);
		pp = p.p + (p.u * dcyl);
		if (pp.z() < DBL_EPSILON)	// столкновение с наружным торцом означает вылет
			return DBL_MAX;
		else
			return dcyl;
	}
	else
	{
		double dprism = mcGeometry::getDistanceToPrismOutside(p.p, p.u, a_, a_, z_) + DBL_EPSILON;
		if (dprism == DBL_MAX)
			return DBL_MAX;

		geomVector3D pp = p.p + (p.u * dprism);

		double xx = pp.x(), yy = pp.y(), zz = pp.z();
		if (zz > z_ / 2 || xx < x1_ || xx > x2_ || yy < y1_ || yy > y2_)
			return dprism;

		// ѕеребираем круги на предмет нахождени€ в каком либо.
		// ≈сли так, то тогда точка пересечени€ будет точкой пересечени€ изнутри цилиндра.
		// ¬ противном случае оп€ть возвращаем dprism.
		// ƒл€ ускорени€ сначала определе€ем индексы пр€моугольника охватывающего точку и образованного центрами отверстий.

		unsigned i, j;
		for (i = 0; i < xs_.size(); i++) if (xs_[i] >= xx) break;
		for (j = 0; j < ys_.size(); j++) if (ys_[j] >= yy) break;

		unsigned ii = 1000, jj = 1000;
		if (i > 0 && j > 0)
		{
			double fx = xx - xs_[i - 1], fy = yy - ys_[j - 1];
			if (fx*fx + fy*fy <= ds_[i - 1] * ds_[i - 1] / 4) { ii = i - 1; jj = j - 1; }
		}
		if (ii == 1000 && j > 0 && i < xs_.size())
		{
			double fx = xx - xs_[i], fy = yy - ys_[j - 1];
			if (fx*fx + fy*fy <= ds_[i] * ds_[i] / 4) { ii = i; jj = j - 1; }
		}
		if (ii == 1000 && i > 0 && j < ys_.size())
		{
			double fx = xx - xs_[i - 1], fy = yy - ys_[j];
			if (fx*fx + fy*fy <= ds_[i - 1] * ds_[i - 1] / 4) { ii = i - 1; jj = j; }
		}
		if (ii == 1000 && i < xs_.size() && j < ys_.size())
		{
			double fx = xx - xs_[i], fy = yy - ys_[j];
			if (fx*fx + fy*fy <= ds_[i] * ds_[i] / 4) { ii = i; jj = j; }
		}

		if (ii != 1000)
		{
			pp(0) -= xs_[ii]; pp(1) -= ys_[jj];
			double dcyl = mcGeometry::getDistanceToCylinderInside(pp, p.u, ds_[ii] / 2, hs_[jj]);
			return dcyl + dprism;
		}

		return dprism;
	}
}

double mcPTLasVegas::getDNearInside(const geomVector3D& p) const
{
	// ѕропускаем ускорение за счет dnear
	return 0;

	//// —начала определ€ем минимальное рассто€ние до параллелепипеда.
	//// «атем перебираем 4 соседних цилиндра.
	//// ¬ идоге выбираем меньшее.
	//double x = p.x(), y = p.y(), z = p.z();
	//double dx = a_ / 2 - fabs(x), dy = a_ / 2 - fabs(y), dz = MIN(z_ - z, z);
	//double dppd = MIN(MIN(dx, dy), dz);

	//unsigned i, j;
	//for (i = 0; i < xs_.size(); i++) if (xs_[i] >= x) break;
	//for (j = 0; j < ys_.size(); j++) if (ys_[j] >= y) break;

	//double d00 = DBL_MAX, d01 = DBL_MAX, d10 = DBL_MAX, d11 = DBL_MAX;
}

void mcPTLasVegas::dump(ostream& os) const
{
	mcTransport::dump(os);
}

void mcPTLasVegas::dumpVRML(ostream& os) const
{
	os << "# Cylindrical grid: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	for (unsigned j = 0; j < ys_.size(); j++)
	{
		for (unsigned i = 0; i < xs_.size(); i++)
		{
			dumpVRMLCylinder(os, ds_[i] / 2, 0, hs_[j], xs_[i], ys_[j]);
		}
	}

	dumpVRMLPrism(os, a_, a_, z_);

	os << "  ]" << endl;
	os << "}" << endl;
}
