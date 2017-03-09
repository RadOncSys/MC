// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// Типы распределений источников
enum mc_distr_t { MCP_CONST_DENS = 0, MCP_GAUSSIAN, MCP_WATERBAG, MCP_K_V, MCP_ARBITRARY };


class mcSourceDistributed : public mcSource
{
public:
	mcSourceDistributed(void);
	mcSourceDistributed(const char* name, int nThreads,
		mc_particle_t type, double ke, const geomVector3D& p, const geomVector3D& v,
		mc_distr_t distr, double Emitx, double rbeamx, double beamAnglex, double Emity,
		double rbeamy, double beamAngley);

	void init(mc_particle_t type       // тип частиц
		, double ke               // кинетическая энергия
		, const geomVector3D& p   // точка рождения частиц
		, const geomVector3D& v   // направление движения (единичный вектор)
		, mc_distr_t distr        // тип распределения
		, double Emitx  // х-эмиттанс, в см. "Нормализованный", без пи, т.е. просто rbeam*rbeta, rbeta = vrmax/C, C - скорость света
		, double rbeamx // полуось по х
		, double beamAnglex // угловое расхождение по х (из-за оптики, не связано с эмиттансом)
		, double Emity            //
		, double rbeamy           //  все то же самое для y
		, double beamAngley       //
	);

	void sample(mcParticle& p, mcThread* thread) override;

	void dumpVRML(ostream& os) const override;

	friend ostream& operator << (ostream& os, const mcSourceDistributed& s)
	{
		os << (const mcSource&)s;
		os << "TYPE = \t" << s.type_ << endl;
		os << "KE = \t" << s.ke_ << endl;
		os << "POSITION = \t" << s.p_ << endl;
		os << "DIRECTION = \t" << s.v_ << endl;
		return os;
	}

protected:
	mc_distr_t distr_;  // тип распределения
	double Emitx_;  // х-эмиттанс, в см. "Нормализованный", без пи, т.е. просто rbeam*rbeta, rbeta = vrmax/C, C - скорость света
	double rbeamx_; // полуось по х
	double rbetax_; // Tmitx/rbeamx
	double beamAnglex_; // угловое расхождение по х (из-за оптики, не связано с эмиттансом)
	double Emity_;            //
	double rbeamy_;           //  все то же самое для y
	double rbetay_;           //
	double beamAngley_;       //
	double beta0_;   // Vz/C, sqrt(ke/Ep/(ke/Ep+1)) - с учетом релятивизма

	mc_particle_t type_;
	double ke_;       // в МэВ                (проверить!!!!!!!!!!!!!!!!!!!)
	geomVector3D p_;  // координаты центра источника
	geomVector3D v_;  // напрвление инжекции, пока допустимо только (0,0,1)
	int q_; // заряд
	// для произвольного распределения - указатель на внешнюю функцию, возвращающую double от 0 до 1
	//  double (*fdistr_)(double x,double y,double px,double py);
};
