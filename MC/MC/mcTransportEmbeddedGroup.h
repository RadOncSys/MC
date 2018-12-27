// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// ����� ���������� � � ������ ���������������� ��������, ��������� � ������ ������.
// ��� �� ���� ������ ������ ���������.
// �� �� �������� ��������, � ������� �������� ��� �������, ���������� �� ��������� ��������.
// ���������������� beginTransport ���������� ���� �� ��������� � ������ ��������� ������
// � ��� ������� �������� ������� � ����� ��������� ������.
// ���� ���, �� �������� ������� ���� ������������ �������, ���� ���������� ��� ������.
class mcTransportEmbeddedGroup : public mcTransport
{
public:
	mcTransportEmbeddedGroup(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportEmbeddedGroup(void);

	void beginTransportInside(mcParticle& p) override;
	void endTransport(mcParticle* particle) override;

	/// <summary>
	/// ���������� ������� ���������� ����������.
	/// </summary>
	void addTransport(mcTransport* t);

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	mcTransport* getInternalTransportByName(const char* name) override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	std::vector<mcTransport*> embeddedTransports_;
};
