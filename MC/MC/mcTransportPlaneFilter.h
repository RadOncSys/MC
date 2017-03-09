// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// Класс специального фильтра.
// Идея - пропускание частиц через бесконечную плоскость 
// (вероятно перпендикулярную оси Z) и сбора информации о них.
// В зависимости от установленных ссылок на предыдущий и последующий транспорты 
// возможно испльзование в качестве поглотителя частиц.
class mcTransportPlaneFilter : public mcTransport
{
public:
	mcTransportPlaneFilter(void);
	mcTransportPlaneFilter(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportPlaneFilter(void);

	// Начало транспорта переписано, так как не нужно запускать симуляцию частиц.
	// Их просто нужно зарегистрировать и пропустить дальше.
	void beginTransport(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;
};
