// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// Класс источника C60 для дистанционных аппаратов.
// Данный класс в отличие от других обеспечивает транспорт внутри самого источника.
// Возможно правильное решение - выделить материалы в класс транспорта и обеспечить
// механизм запуск стартуя изнутри объекта.
//
// Автор: GG 
// Дата: 09.02.2013
//-----------------------------------------------------------------------------------------
class mcSourceCylindricalC60 : public mcSource
{
public:
	mcSourceCylindricalC60(const char* name, int nThreads,
		const geomVector3D& p, const geomVector3D& v, double r, double h);
	virtual ~mcSourceCylindricalC60(void);

	void sample(mcParticle& p, mcThread* thread) override;

	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceCylindricalC60& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \tC60" << endl;
		os << "NAME = \t" << s.getName() << endl;
		os << "DIAMETER = \t" << 2 * s.r_ << endl;
		os << "HEIGHT = \t" << s.h_ << endl;
		os << "POSITION = \t" << s.p_ << endl;
		return os;
	}

protected:
	geomVector3D p_;
	geomVector3D v_;
	double r_;
	double h_;
};
