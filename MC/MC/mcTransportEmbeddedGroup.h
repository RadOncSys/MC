// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс транспорта в в группе непересекающихся объектов, вложенных в данный объект.
// Сам по себе данный объект прозрачен.
// Но он является объектом, в который попадают все частицы, вылетающие из вложенных объектов.
// Переопределенный beginTransport определяет есть ли попадание в другой вложенный объект
// и при наличии передает частицу в такой ближайший объект.
// если нет, то передает частицу либо охватвающему объекту, либо следующему как обычно.
class mcTransportEmbeddedGroup : public mcTransport
{
public:
	mcTransportEmbeddedGroup(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportEmbeddedGroup(void);

	void beginTransportInside(mcParticle& p) override;
	void endTransport(mcParticle* particle) override;

	/// <summary>
	/// Добавление объекта вложенного транспорта.
	/// </summary>
	void addTransport(mcTransport* t);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	mcTransport* getInternalTransportByName(const char* name) override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	std::vector<mcTransport*> embeddedTransports_;
};
