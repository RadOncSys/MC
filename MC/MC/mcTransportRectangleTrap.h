// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта который поглощает все частицы, за исключением тех
// которые летят в положительном направлении и попадают в заданный прямоугольник.
class mcTransportRectangleTrap : public mcTransport
{
public:
	mcTransportRectangleTrap(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	~mcTransportRectangleTrap(void);

	// Начало транспорта переписано, так как не нужно запускать симуляцию частиц.
	// Их просто нужно отфильтровать.
	void beginTransport(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;

	void SetFieldSize(double x1, double x2, double y1, double y2);

protected:
	double fsx1_;
	double fsx2_;
	double fsy1_;
	double fsy2_;
};
