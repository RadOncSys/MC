#include ".\mcgeomtriangleside.h"
#include <float.h>

mcGeomTriangleSide::
mcGeomTriangleSide(const geomVector3D& p
	, const geomVector3D& Vx
	, const geomVector3D& Vy
	, double ax, double ay)
	:P_(p)
	, Vx_(Vx)
	, Vy_(Vy)
	, ax_(ax)
	, ay_(ay)
{
	Vx_.normalize();
	Vy_.normalize();
	N_ = Vx_^Vy_;
	N_.normalize();
	VRy_ = N_^Vx_;

	// ������ ������������
	p_[0].set(0, 0);
	p_[1].set(ax_, 0);
	p_[2].set(ay*(Vx_*Vy_), ay*(VRy_*Vy_));

	// ������� � ������� �����
	for (int i = 0; i < 3; i++)
	{
		v_[i] = p_[(i + 1) % 3] - p_[i];
		v_[i].normalize();
		n_[i] = v_[i];
		n_[i].turnLeft();
	}

	// ����� �����
	a_[0] = ax_;
	a_[1] = ((Vx_*ax_) - (Vy_*ay_)).length();
	a_[2] = ay_;
}

double mcGeomTriangleSide::getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const
{
	// ������ ����� ��� ����������, � ������� ��������� �������������
	double h = (P_ - p)*N_;

	// ������ ������� ����������� � ���������
	double cosn = v*N_;
	if (fabs(cosn) < DBL_EPSILON)
		return DBL_MAX; // ����������� ���������

	if (fabs(h) < DBL_EPSILON && ((inside && cosn < 0) || (!inside && cosn > 0)))
		return DBL_MAX; // �������� �� �����������

	// ����� ����������� ������� ����������� � ����������
	double distance = h / cosn;
	if (distance < 0)
		return DBL_MAX; // �������� �� ���������
	geomVector3D pc = p + (v*distance);

	// ���������� ����� �������� � ������� ������������
	geomVector3D pr = pc - P_;
	geomVector2D R(pr*Vx_, pr*VRy_);

	// ������ � ������
	int i;
	double hr[3];
	bool isInside = true;
	for (i = 0; i < 3; i++)
	{
		hr[i] = (R - p_[i]) * n_[i];
		if (hr[i] < 0)
			isInside = false;
	}
	if (isInside)
		return distance;
	else
		return DBL_MAX; // � ������ ����������� �� ���� ������������
}

double mcGeomTriangleSide::getDNear(const geomVector3D& p) const
{
	// ������ ����� ��� ����������, � ������� ��������� �������������
	double h = (p - P_)*N_;

	// �������� ����� �� ��������� �������������� (�� �������)
	geomVector3D pp = p + (N_*h);

	// ���������� ����� �������� � ������� ������������
	geomVector3D pr = pp - P_;
	geomVector2D R(pr*Vx_, pr*VRy_);

	// ������ � ������
	int i;
	double hr[3];
	bool isInside = true;
	for (i = 0; i < 3; i++)
	{
		hr[i] = (R - p_[i]) * n_[i];
		if (hr[i] < 0)
			isInside = false;
	}
	// ���� �������� � ����� ��������������, �� ����������� ����������� ����� ������
	if (isInside)
		return fabs(h);

	// ����������� ����� ����������� ���������� �� ������ �� �����.
	// ���� ���������� ����� �� ������� ���������� ������,
	// �� ������ ����������� ����� ����� ����� �����������
	// ����������� �� ������������.
	// � ��������� ������ ����������� ����� ���������� �� ����� �� ������.
	double dmin = DBL_MAX;
	for (i = 0; i < 3; i++)
	{
		double l = (R - p_[i]) * v_[i];
		if (l >= 0 && l <= a_[i] && dmin > hr[i])
			dmin = hr[i];
	}
	if (dmin < DBL_MAX)
		return sqrt(dmin*dmin + h*h);

	for (i = 0; i < 3; i++)
	{
		double d = (R - p_[i]).length();
		if (d < dmin)
			dmin = d;
	}
	return sqrt(dmin*dmin + h*h);
}
