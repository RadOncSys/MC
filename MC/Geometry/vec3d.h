// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mtrx3d.h"

class geomVector3D
{
public:
	geomVector3D() { set(0, 0, 0); }
	geomVector3D(double x, double y, double z) { set(x, y, z); }

	double x() const { return p_[0]; }
	double y() const { return p_[1]; }
	double z() const { return p_[2]; }

	void set(double x, double y, double z) {
		p_[0] = x;	p_[1] = y;	p_[2] = z;	p_[3] = 1;
	}

	double longSide() const {
		double lx = fabs(p_[0]), ly = fabs(p_[1]), lz = fabs(p_[2]);
		if (lx < ly) lx = ly;
		if (lx < lz) lx = lz;
		return lx;
	}

	double length() const { return sqrt(sqLength()); }

	double sqLength() const {
		return p_[0] * p_[0] + p_[1] * p_[1] + p_[2] * p_[2];
	}

	double lengthXY() const { return sqrt(sqLengthXY()); }

	double sqLengthXY() const {
		return p_[0] * p_[0] + p_[1] * p_[1];
	}

	void normalize() {
		double f = length();
		p_[0] /= f; p_[1] /= f; p_[2] /= f;
	}

	geomVector3D transformDirection(const geomMatrix3D& m) const {
		geomVector3D v;
		for (int i = 0; i < 3; i++) {
			v(i) = 0;
			for (int j = 0; j < 3; j++)
				v(i) += p_[j] * m.m_[j][i];
		}
		return v;
	}

	bool isInsideBBox(const geomVector3D& b1, const geomVector3D& b2)const;

	double& operator()(int n) { return p_[n]; }
	double  operator()(int n) const { return p_[n]; }

	geomVector3D operator*(double f) const {
		return geomVector3D(p_[0] * f, p_[1] * f, p_[2] * f);
	}

	geomVector3D operator^(const geomVector3D& v) const {
		return geomVector3D(p_[1] * v.z() - p_[2] * v.y()
			, p_[2] * v.x() - p_[0] * v.z()
			, p_[0] * v.y() - p_[1] * v.x());
	}

	geomVector3D operator+(const geomVector3D& v) const {
		return geomVector3D(p_[0] + v.x(), p_[1] + v.y(), p_[2] + v.z());
	}

	geomVector3D operator-(const geomVector3D& v) const {
		return geomVector3D(p_[0] - v.x(), p_[1] - v.y(), p_[2] - v.z());
	}

	double operator*(const geomVector3D& v) const {
		return p_[0] * v.x() + p_[1] * v.y() + p_[2] * v.z();
	}

	geomVector3D operator*(const geomMatrix3D& m) const {
		geomVector3D v;
		for (int i = 0; i < 4; i++) {
			v(i) = 0;
			for (int j = 0; j < 4; j++)
				v(i) += p_[j] * m.m_[j][i];
		}
		return v;
	}

	void operator*=(double f) {
		p_[0] *= f; p_[1] *= f; p_[2] *= f;
	}

	void operator/=(double f) {
		p_[0] /= f; p_[1] /= f; p_[2] /= f;
	}

	void operator+=(const geomVector3D& v) {
		p_[0] += v.x(); p_[1] += v.y(); p_[2] += v.z();
	}

	void operator-=(const geomVector3D& v) {
		p_[0] -= v.x(); p_[1] -= v.y(); p_[2] -= v.z();
	}

	bool operator==(const geomVector3D& v) const {
		return p_[0] == v.x() && p_[1] == v.y() && p_[2] == v.z();
	}

	inline void operator=(const geomVector3D& v) {
		p_[0] = v.x(); p_[1] = v.y(); p_[2] = v.z();
	}

	friend istream& operator >> (istream&, geomVector3D&);
	friend ostream& operator << (ostream&, const geomVector3D&);

	double p_[4];
};
