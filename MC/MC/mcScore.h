// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcParticle.h"
#include "mcObj.h"

class mcScore : public mcObj
{
public:
	mcScore(const char* module_name, int nThreads);
	virtual ~mcScore(void);

	virtual void ScoreFluence(const mcParticle& particle) { }

	virtual void ScorePoint(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0) {
		etotal_[iThread] += edep;
	}

	virtual void ScoreLine(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0
		, const geomVector3D& p1) {
		etotal_[iThread] += edep;
	}

	// Преобразование зачений энерги в значения дозы
	virtual void CE2D() { }
	bool isConvertedToDose() const { return dconverted_; }

	const char* getModuleName() const { return attached_module_; }

	void setTransport(mcTransport* t) { transport_ = t; }

	// Вывод изображения детектора в формате VRML.
	// Матрица преобразования координат в мировую систему необходима,
	// так как scoring задан в системе координат объекта транспорта, с которым он связан.
	virtual void dumpVRML(ostream&) const;

	// Вывод содержимого перемещен не в оператор, а в функцию,
	// так как требуется поддержка наследования.
	virtual void dumpStatistic(ostream&) const;

	void setDensity(double d) { density_ = d; }

	double etotal() const;

protected:
	char attached_module_[64];
	int nThreads_;
	bool dconverted_;
	double* etotal_;
	double density_;
	mcTransport* transport_;
};
