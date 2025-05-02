// Radiation Oncology Monte Carlo open source project
//
// Author: [2025] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"
#include <vector>

// Источник в виде радиально симметричного (вокруг оси Z) потока фотонов.
// Частицы беруттся коллекции предыдущих симуляций циклически и поворачиваются на случайный угол.
// Частицы появляются на поверхности того объема, который отвечал за их сбор при предыдущих симуляциях.
// Этот объем в норме небольшой и использующая программа должна обеспечивать 
// нахождение этого объема полностью внутри содержащего источник объекта, 
// например, воздушной полостии эндостата в брахитерапии.
class mcEBTSourcePhsp : public mcSource
{
public:
	mcEBTSourcePhsp(const char* name, int nThreads, double z);

	void sample(mcParticle& p, mcThread* thread) override;
	//void loadData(istream& is);
	void readFromMemory(void* buffer);
	void dumpVRML(ostream& os) const override;

protected:
	// Смещение источника (содержащего частицы в собственной системе координат) по оси Z
	double z_;
	mc_particle_t t_;
	int nparticles_;
	int npars_;
	std::vector<int> idx_;
	std::vector<double> phsp_;
};
