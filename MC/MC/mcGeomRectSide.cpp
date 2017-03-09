#include ".\mcgeomrectside.h"
#include <float.h>

mcGeomRectSide::
mcGeomRectSide(const geomVector3D& p
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
}

double mcGeomRectSide::getDistance(const geomVector3D& p, const geomVector3D& v, bool inside) const
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

	// ���������� ����� �������� � ������� ��������������
	geomVector3D pr = pc - P_;
	double x = pr*Vx_;
	double y = pr*Vy_;

	// ���� ����� ����������� �� � �����, �� ���������� �����������
	if (x >= 0 && x <= ax_ && y >= 0 && y <= ay_)
		return distance;
	else
		return DBL_MAX; // � ������ ����������� �� ���� ��������������
}

double mcGeomRectSide::getDNear(const geomVector3D& p) const
{
	// ������ ����� ��� ����������, � ������� ��������� �������������
	double h = (p - P_)*N_;

	// �������� ����� �� ��������� �������������� (�� �������)
	geomVector3D pp = p + (N_*h);

	// ���������� ����� �������� � ������� ��������������
	geomVector3D pr = pp - P_;
	double x = pr*Vx_;
	double y = pr*Vy_;

	// ���� �������� � ����� ��������������, �� ����������� ����������� ����� ������
	if (x >= 0 && x <= ax_ && y >= 0 && y <= ay_)
		return fabs(h);

	// ������� ��������
	if (x < 0)
	{
		if (y < 0) return sqrt(x*x + y*y + h*h);
		else if (y < ay_) return sqrt(x*x + h*h);
		else return sqrt(x*x + (y - ay_)*(y - ay_) + h*h);
	}
	else if (x > ax_)
	{
		if (y < 0) return sqrt((x - ax_)*(x - ax_) + y*y + h*h);
		else if (y < ay_) return sqrt((x - ax_)*(x - ax_) + h*h);
		else return sqrt((x - ax_)*(x - ax_) + (y - ay_)*(y - ay_) + h*h);
	}
	else
	{
		if (y < 0) return sqrt(y*y + h*h);
		else return sqrt((y - ay_)*(y - ay_) + h*h);
	}
}
