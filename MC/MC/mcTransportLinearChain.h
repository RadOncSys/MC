// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

//  ласс оболочки над линейной цепочкой. 
// ¬нутренние модули обмениваютс€ частицами в обычном отсортированном по Z стиле.
// ѕри окончании транспорта первым или последним модулем транспорта цепочки
// „астица передаетс€ данному классу оболочки подобно тому как работают вложенные модули.
// ѕри получении частицы дл€ транспорта сам отыскивает модкль цепочки, куда частица попадет в первую очередь.
class mcTransportLinearChain : public mcTransport
{
public:
	mcTransportLinearChain(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportLinearChain(void);

	void beginTransport(mcParticle& p) override;
	void beginTransportInside(mcParticle& p) override;

	double getDistanceOutside(mcParticle& p) const override;

	mcTransport* getInternalTransportByName(const char* name) override;

	// ƒобавление объекта транспорта цепочки.
	// –егистраци€ участников транспорта необходима дл€ определени€ 
	// в какой транспорт попадет частица, когда она проникнет в цепочку снаружи.
	void addTransport(mcTransport* t);

	// ¬ызываетс€ по окончании парсинга дл€ заверешени€ установки всех св€зей
	void completeInit();

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	std::vector<mcTransport*> chainTransports_;
};
