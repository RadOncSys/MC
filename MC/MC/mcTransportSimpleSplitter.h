// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// ����� ���������� ���� ������������ �������.
// � ������� ����� ������ �� ���������� ��������������.
// ������ �������, ������������ ��������� � ������������� ����������� ������������
// �� ��������� ����� � ��������������� ������������ ������.
class mcTransportSimpleSplitter : public mcTransport
{
public:
	mcTransportSimpleSplitter(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, mc_particle_t ptype, int nsplit);
	~mcTransportSimpleSplitter(void);

	// ������ ���������� ����������, ��� ��� �� ����� ��������� ��������� ������.
	// �� ������ ����� ���������������� � ���������� ������.
	void beginTransport(mcParticle& p) override;
	void beginTransportInside(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;

protected:
	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	enum mc_particle_t ptype_;
	int nsplit_;
};
