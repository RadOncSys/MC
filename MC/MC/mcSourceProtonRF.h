// Radiation Oncology Monte Carlo open source project
//
// Author: [2025] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"
#include <vector>

// »сточник в виде набора протонов на входе в оптическую скамью.
// ќбычно набор представл€етс€ текстовой таблицей, расчитанной в радиофизических симул€ци€х.
// ‘айл частиц должен удовлетвор€ть простому текстовому формату:
// X[mm] Vx[rad] Y[mm] Vy[rad] dE/E[%]
class mcSourceProtonRF : public mcSource
{
public:
	mcSourceProtonRF(const char* name, int nThreads, double z);

	void sample(mcParticle& p, mcThread* thread) override;
	void loadData(istream& is);
	void dumpVRML(ostream& os) const override;

protected:
	// —мещение источника (содержащего частицы в собственной системе координат) по оси Z
	double z_;
	double e0_;
	int nparticles_;
	std::vector<int> idx_;
	std::vector<mcRng> rng_;
	std::vector<mcParticle> particles_;
};
