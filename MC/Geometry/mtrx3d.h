// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <math.h>
#include <iostream>
using namespace std;

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

class geomVector3D;

class geomMatrix3D
{
public:
	geomMatrix3D();
	geomMatrix3D(double a00, double a01, double a02, double a03,
		double a10, double a11, double a12, double a13,
		double a20, double a21, double a22, double a23,
		double a30, double a31, double a32, double a33);

	geomMatrix3D(const geomMatrix3D& matrix) {
		memcpy(m_, matrix.m_, 16 * sizeof(double));
	}

	void makeInverse();
	void makeUnit();

	// Матрицы для получения преобразований координат:
	// 1) Параллельный перенос.
	static geomMatrix3D ParallelShift(double Ax, double Ay, double Az) {
		return geomMatrix3D(1., 0., 0., 0.,
			0., 1., 0., 0.,
			0., 0., 1., 0.,
			Ax, Ay, Az, 1.);
	}
	// 2) Масштабирование.
	static geomMatrix3D Scaling(double S) {
		return geomMatrix3D(S, 0., 0., 0.,
			0., S, 0., 0.,
			0., 0., S, 0.,
			0., 0., 0., 1.);
	}
	static geomMatrix3D Scaling(double Sx, double Sy, double Sz) {
		return geomMatrix3D(Sx, 0., 0., 0.,
			0., Sy, 0., 0.,
			0., 0., Sz, 0.,
			0., 0., 0., 1.);
	}
	// 3) Вращение вокруг оси Ох на угол fi1:
	static geomMatrix3D RotationAroundOx(double fi1) {    // fi1 в радианах
		return geomMatrix3D(1., 0., 0., 0.,
			0., cos(fi1), sin(fi1), 0.,
			0., -sin(fi1), cos(fi1), 0.,
			0., 0., 0., 1.);
	}
	static geomMatrix3D RotationAroundOx(int fi1) {       // fi1 в градусах
		double f1 = fi1*PI / 180.;
		return geomMatrix3D(1., 0., 0., 0.,
			0., cos(f1), sin(f1), 0.,
			0., -sin(f1), cos(f1), 0.,
			0., 0., 0., 1.);
	}
	// 4) Вращение вокруг оси Оy на угол fi2:
	static geomMatrix3D RotationAroundOy(double fi2) {    // fi2 в радианах
		return geomMatrix3D(cos(fi2), 0., -sin(fi2), 0.,
			0., 1., 0., 0.,
			sin(fi2), 0., cos(fi2), 0.,
			0., 0., 0., 1.);
	}
	static geomMatrix3D RotationAroundOy(int fi2) {       // fi2 в градусах
		double f2 = fi2*PI / 180.;
		return geomMatrix3D(cos(f2), 0., -sin(f2), 0.,
			0., 1., 0., 0.,
			sin(f2), 0., cos(f2), 0.,
			0., 0., 0., 1.);
	}
	// 5) Вращение вокруг оси Оz на угол fi3:
	static geomMatrix3D RotationAroundOz(double fi3) {    // fi3 в радианах
		return geomMatrix3D(cos(fi3), sin(fi3), 0., 0.,
			-sin(fi3), cos(fi3), 0., 0.,
			0., 0., 1., 0.,
			0., 0., 0., 1.);
	}
	static geomMatrix3D RotationAroundOz(int fi3) {       // fi3 в градусах
		double f3 = fi3*PI / 180.;
		return geomMatrix3D(cos(f3), sin(f3), 0., 0.,
			-sin(f3), cos(f3), 0., 0.,
			0., 0., 1., 0.,
			0., 0., 0., 1.);
	}
	// 6) Матрица вращения из однрй системы координат в другую.
	//    Аргументами являются единичные векторы осей новой 
	//    системы в старой системе координат
	static geomMatrix3D BuildFromAxis(const geomVector3D& ax
		, const geomVector3D& ay
		, const geomVector3D& az);

	double& operator()(int m, int n) {
		return m_[m][n];
	}

	geomMatrix3D operator*(const geomMatrix3D& m) const {
		geomMatrix3D	rmtrx;
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int k = 0; k < 4; k++)
					rmtrx.m_[i][j] += m_[i][k] * m.m_[k][j];
		return rmtrx;
	}

	bool operator==(const geomMatrix3D& m) const {
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				if (m_[i][j] != m.m_[i][j])
					return false;
		return true;
	}

	friend istream& operator >> (istream&, geomMatrix3D&);
	friend ostream& operator << (ostream&, const geomMatrix3D&);

	double m_[4][4];
};
