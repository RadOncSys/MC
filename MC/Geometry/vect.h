// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <math.h>

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

class CVect_int
{
public:
	CVect_int() :x(0), y(0) {}

	CVect_int(int abs, int ord) { x = abs; y = ord; }

	CVect_int set(int abs, int ord) {
		x = abs; y = ord;
		return *this;
	}

	// Поворот против часовой стрелке на угол 45·а	в физических координатах
	// и по часовой стрелке - в экранных
	void turn(int a) {
		int abs[8] = { 0,  1,1,1,0,-1,-1,-1 };
		int ord[8] = { -1,-1,0,1,1, 1, 0, -1 };
		int i;
		for (i = 0; i < 8; i++) if (x == abs[i] && y == ord[i]) break;
		i = i + a;
		if (i > 7) i = i - 8;
		if (i < 0) i = i + 8;
		x = abs[i];
		y = ord[i];
	}

	// Превращает вектор в почти так же ориентированный ветор
	// с координатами не больше 1
	void unitv()
	{
		double abs, ord;
		abs = (x < 0) ? -x : x;
		ord = (y < 0) ? -y : y;
		if (abs == 0) {
			x = 0;
			y = (y > 0) ? 1 : -1;
			return;
		}
		if (ord / abs < 0.58) {
			y = 0;
			x = (x > 0) ? 1 : -1;
			return;
		}
		if (ord / abs > 1.73) {
			x = 0;
			y = (y > 0) ? 1 : -1;
			return;
		}
		x = (x > 0) ? 1 : -1;
		y = (y > 0) ? 1 : -1;
		return;
	}

	// Возвращает расстояние между this вектором и вектором аргументом
	double dist(const CVect_int& one)const {
		return sqrt(double((one.x - x)*(one.x - x) + (one.y - y)*(one.y - y)));
	}

	int sqDist(const CVect_int& one)const {
		return (one.x - x)*(one.x - x) + (one.y - y)*(one.y - y);
	}

	// Возвращает расстояние от this вектора до отрезка Р1-Р2
	double dist(const CVect_int& one, const CVect_int& two)const {
		double l1, l2, l3, h;
		l1 = one.dist(two);
		l2 = this->dist(two);
		l3 = this->dist(one);
		h = (l1*l1 + l3*l3 - l2*l2) / l1 / 2;
		return sqrt(l3*l3 - h*h);
	}

	// возвращает длину вектора
	double length() const
	{
		return sqrt(double(x*x + y*y));
	}

	// Определяет, находится ли this вектор внутри прямоугольника
	bool inside(int left, int top, int right, int bottom)const {
		if (top < bottom) { int t = bottom; bottom = top; top = t; }
		return (x >= left && x <= right && y <= top && y >= bottom);
	}

	// Если this вектор вне прямоугольника RECT перемещает вектор
	// в ближайшую точку на границе прямоугольника
	void inside_move(int left, int top, int right, int bottom) {
		if (top < bottom) { int t = bottom; bottom = top; top = t; }
		x = (x < left) ? left : x;
		x = (x > right) ? right : x;
		y = (y < bottom) ? bottom : y;
		y = (y > top) ? top : y;
	}

	// Возвращает угол от -pi до +pi
	double V_angle() const
	{
		double b;
		b = sqrt(double(x*x + y*y));
		if (x == 0) return (y > 0) ? PI / 2 : -PI / 2;
		if (y == 0) return (x > 0) ? 0 : PI;
		if (x > 0) return asin(y / b);
		return (y > 0) ? PI - asin(y / b) : -PI - asin(y / b);
	}

	// Превращает double координаты в ближайший целочисленный узел
	void to_int(double xd, double yd) {
		int _x;
		int _y;
		_x = (int)xd;
		_y = (int)yd;
		_y = (yd - _y >= 0.5) ? _y + 1 : _y;
		_x = (xd - _x >= 0.5) ? _x + 1 : _x;
		x = _x;
		y = _y;
	}

	CVect_int operator* (int k) const
	{
		return CVect_int(k*x, k*y);
	}

	CVect_int operator+ (const CVect_int& one) const
	{
		return CVect_int(x + one.x, y + one.y);
	}

	bool operator== (const CVect_int& one) const {
		return x == one.x&&y == one.y;
	}

	bool operator!=(const CVect_int& one) const {
		return !(x == one.x&&y == one.y);
	}

	CVect_int operator=(const CVect_int& one) {
		x = one.x;
		y = one.y;
		return *this;
	}

	CVect_int operator-(const CVect_int& one) const
	{
		CVect_int a(x - one.x, y - one.y);
		return a;
	}

	CVect_int operator-() const
	{
		CVect_int a(-x, -y);
		return a;
	}

public:
	int x, y;
};
