// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"


// TODO: Класс стартовал с существующего транспорта mcPTLasVegas
// похожего своей матрицей отверстией. 
// Теперь его функционал нужно переписать под гребенчатый фильтр.



// Класс транспорта в гребенчатом фильтре, представляущем 2D матрицу пирамид, 
// состоящих из корпичиков квадратного сечения одинаковой высоты, но разного размера.
// Изначально класс разрабатывался для транспорта протонов в проекте ОКО.
class mcTransportGridFilter : public mcTransport
{
public:
	mcTransportGridFilter(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportGridFilter(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	// Нужно для тестирования
	double A() const { return a_; }
	double Z() const { return z_; }
	const std::vector<double>& Ds() const { return ds_; }
	const std::vector<double>& Hs() const { return hs_; }
	const std::vector<double>& Xs() const { return xs_; }
	const std::vector<double>& Ys() const { return ys_; }

protected:
	// Объект описывается как двумерная сетка цилиндрических отверстий (5х5)
	double a_;	// сторона квадратного фантома
	double z_;	// толщина фантома

	// Диаметры отверстий
	std::vector<double> ds_;

	// Глубины отверстий
	std::vector<double> hs_;

	// Позиции узлов сетки отверстий
	std::vector<double> xs_;
	std::vector<double> ys_;

	// Служебные переменные, определяющие где могут находиться высверленные отверстия
	double x1_, x2_, y1_, y2_, z2_;
	double ax_, ay_;	// симметричный прямоугольник, охватывающий область отверстий
};
