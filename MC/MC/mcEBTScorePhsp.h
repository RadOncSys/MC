// Radiation Oncology Monte Carlo open source project
//
// Author: [2025] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс регистрации фотонов от рентгеновской трубки электронной брахитерапии. 
// Продуктом работы является текстовый список частиц в пространстве 
// в норме на поверхности объекта, к которому прикреплен данный скоринг. 
// Частицы регистрируются как есть, т.е. нет предположений типа радиальной симметрии.
// Симметрия реализуется при необходимости на уровне использующего результаты скоринга 
// источника элементарным поворотом вокруг оси симметрии на случайный угол.

class mcEBTScorePhsp : public mcScore
{
public:
	mcEBTScorePhsp(const char* module_name, int nThreads, const char* outfname);
	virtual ~mcEBTScorePhsp();

	void ScoreFluence(const mcParticle& particle) override;

	void SavePHSP() const;

	void dumpVRML(ostream&) const override;
	void dumpStatistic(ostream&) const override;

protected:
	int npars_;
	std::string outfname_;
	std::vector<std::vector<double>> phsp_;
};
