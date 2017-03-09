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
	mcTransportRectangleTrap(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double width, double length);
	~mcTransportRectangleTrap(void);

	// Начало транспорта переписано, так как не нужно запускать симуляцию частиц.
	// Их просто нужно отфильтровать.
	void beginTransport(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;

protected:
	double swidth_;
	double slength_;
};
