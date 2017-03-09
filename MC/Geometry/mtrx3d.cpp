#include "mtrx3d.h"
#include "vec3d.h"
#include "text.h"

geomMatrix3D::
geomMatrix3D(double a00, double a01, double a02, double a03,
	double a10, double a11, double a12, double a13,
	double a20, double a21, double a22, double a23,
	double a30, double a31, double a32, double a33)
{
	m_[0][0] = a00;	m_[0][1] = a01;	m_[0][2] = a02;	m_[0][3] = a03;
	m_[1][0] = a10;	m_[1][1] = a11;	m_[1][2] = a12;	m_[1][3] = a13;
	m_[2][0] = a20;	m_[2][1] = a21;	m_[2][2] = a22;	m_[2][3] = a23;
	m_[3][0] = a30;	m_[3][1] = a31;	m_[3][2] = a32;	m_[3][3] = a33;
}

geomMatrix3D::geomMatrix3D()
{
	memset(m_, 0, 16 * sizeof(double));
}

void geomMatrix3D::makeInverse()
{
	double r1[8], r2[8], r3[8], r4[8];
	double *s[4], *tmprow;

	s[0] = &r1[0];
	s[1] = &r2[0];
	s[2] = &r3[0];
	s[3] = &r4[0];

	register int i, j, p, jj;
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			s[i][j] = m_[i][j];
			if (i == j) s[i][j + 4] = 1.0;
			else     s[i][j + 4] = 0.0;
		}
	}
	double scp[4];
	for (i = 0; i < 4; i++)
	{
		scp[i] = fabs(s[i][0]);
		for (j = 1; j < 4; j++)
			if (fabs(s[i][j]) > scp[i]) scp[i] = fabs(s[i][j]);
		if (scp[i] == 0.0) return; // singular matrix!
	}

	int pivot_to;
	double scp_max;
	for (i = 0; i < 4; i++)
	{
		// select pivot row
		pivot_to = i;
		scp_max = fabs(s[i][i] / scp[i]);
		// find out which row should be on top
		for (p = i + 1; p < 4; p++)
			if (fabs(s[p][i] / scp[p]) > scp_max)
			{
				scp_max = fabs(s[p][i] / scp[p]); pivot_to = p;
			}
		// Pivot if necessary
		if (pivot_to != i)
		{
			tmprow = s[i];
			s[i] = s[pivot_to];
			s[pivot_to] = tmprow;
			double tmpscp;
			tmpscp = scp[i];
			scp[i] = scp[pivot_to];
			scp[pivot_to] = tmpscp;
		}

		double mji;
		// perform gaussian elimination
		for (j = i + 1; j < 4; j++)
		{
			mji = s[j][i] / s[i][i];
			s[j][i] = 0.0;
			for (jj = i + 1; jj < 8; jj++)
				s[j][jj] -= mji*s[i][jj];
		}
	}
	if (s[3][3] == 0.0) return; // singular matrix!

	//
	// Now we have an upper triangular matrix.
	//
	//  x x x x | y y y y
	//  0 x x x | y y y y 
	//  0 0 x x | y y y y
	//  0 0 0 x | y y y y
	//
	//  we'll back substitute to get the inverse
	//
	//  1 0 0 0 | z z z z
	//  0 1 0 0 | z z z z
	//  0 0 1 0 | z z z z
	//  0 0 0 1 | z z z z 
	//

	double mij;
	for (i = 3; i > 0; i--)
	{
		for (j = i - 1; j > -1; j--)
		{
			mij = s[j][i] / s[i][i];
			for (jj = j + 1; jj < 8; jj++)
				s[j][jj] -= mij*s[i][jj];
		}
	}

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			m_[i][j] = s[i][j + 4] / s[i][i];
}

void geomMatrix3D::makeUnit()
{
	memset(m_, 0, 16 * sizeof(double));
	m_[0][0] = m_[1][1] = m_[2][2] = m_[3][3] = 1;
}

geomMatrix3D geomMatrix3D::BuildFromAxis(const geomVector3D& ax
	, const geomVector3D& ay
	, const geomVector3D& az)
{
	return geomMatrix3D(ax.x(), ay.x(), az.x(), 0,
		ax.y(), ay.y(), az.y(), 0,
		ax.z(), ay.z(), az.z(), 0,
		0, 0, 0, 1);
}

istream& operator >> (istream& is, geomMatrix3D& m)
{
	string line;
	double p[4];
	for (unsigned j = 0; j < 4; j++) {
		getline(is, line, '\n');
		GetFloatArray(line, p, 4);
		for (unsigned i = 0; i < 4; i++)
			m.m_[j][i] = p[i];
	}
	return is;
}

ostream& operator << (ostream& os, const geomMatrix3D& m)
{
	for (unsigned j = 0; j < 4; j++) for (unsigned i = 0; i < 4; i++) { os << m.m_[j][i]; if (i < 3) os << '\t'; else os << endl; }
	return os;
}
