// Radiation Oncology Monte Carlo open source project
//
// Author: [2025] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"
#include "mcGeomRectSide.h"

// Класс отверстия в цилиндре, представленного полигоном  в сечении и вертикальными стенками.
// За исключением формы горизонтального сечения похож на класс mcTransportConicalHole
// Поэтому, в типичной конфигурации расчетов ось Z объекта направлена вдоль -Z.
class mcTransportMantleBlock : public mcTransport
{
public:
	mcTransportMantleBlock(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, 
		double r1, double h, std::vector<double> plgnX, std::vector<double> plgnY);
	virtual ~mcTransportMantleBlock(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

protected:
	bool isPointInPlgn(double x, double y) const;

	double r1_; // внешний радиус
	double h_;  // высота цилиндра
	int nsides_;
	std::vector<double> x_;
	std::vector<double> y_;
	std::unique_ptr<std::vector<mcGeomRectSide>> sides_;
};
