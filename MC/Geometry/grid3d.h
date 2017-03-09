// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "rect3d.h"

class geomGrid3D : public geomRect3D
{
public:
	geomGrid3D();
	geomGrid3D(const geomGrid3D&);

	unsigned
		numCol() const { return nx_; }
	unsigned
		numRow() const { return ny_; }
	double
		pixelSX() const { return psx_; }
	double
		pixelSY() const { return psy_; }
	geomVector3D
		getPointByIndex(int, int)const;
	void
		getPointByIndex(int i, int j, double& x, double& y)const;

	void
		setGrid(int nx, int ny, double psx, double psy);
	void
		setResolution(double x, double y);
	// Интелектуальная установка сетки при которой обязательно присутствуют
	// точки на центральных осях и при этом шаг точно равен указанному.
	// Ячейки сетки, при этом, заполняют всю рамку.
	// Следует обратить внимание, что сетка задает границы ячеек,
	// в то время как дозиметрический сервер считает дозу в центрах ячеек.
	// Например, превая точка расчета определяется ка:
	//   xfirst=m.minX()+0.5*psx
	void
		setAdjustedResolution(double x, double y);
	void
		setSize(int nx, int ny);

	// Индексы по заданным физическим координатам,
	// проверки на границы матрицы нет, возможны отрицательные значения,
	// пользоатель должен сам анализировать результат
	int
		getColIndex(double x)const;
	int
		getRowIndex(double y)const;
	double
		getCol(unsigned i)const;
	double
		getRow(unsigned j)const;

	void scale(double f);

	void operator=(const geomGrid3D&);

	friend istream& operator >> (istream&, geomGrid3D&);
	friend ostream& operator << (ostream&, const geomGrid3D&);

protected:
	unsigned nx_, ny_;		// numbers of cels
	double psx_, psy_; 	// grid cell size
};
