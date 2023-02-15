#include "mcCSNuclear.h"

void mcCSNuclearForAngleSpectrum::Clear()
{
	SourceBinUps.clear();
	SourceSpectrum.clear();
}

void mcCSNuclearParticleEnergy::Clear()
{
	Energy = 0;
	TotalCrossSection = 0;
	ProtonCrossSection = 0;
	NeutronCrossSection = 0;

	ProtonAngles.clear();
	NeutronAngles.clear();
}

void mcCSNuclear::Clear()
{
	Energies.clear();
}
