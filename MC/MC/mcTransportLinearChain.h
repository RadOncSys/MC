// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// ����� �������� ��� �������� ��������. 
// ���������� ������ ������������ ��������� � ������� ��������������� �� Z �����.
// ��� ��������� ���������� ������ ��� ��������� ������� ���������� �������
// ������� ���������� ������� ������ �������� ������� ���� ��� �������� ��������� ������.
// ��� ��������� ������� ��� ���������� ��� ���������� ������ �������, ���� ������� ������� � ������ �������.
class mcTransportLinearChain : public mcTransport
{
public:
	mcTransportLinearChain(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportLinearChain(void);

	void beginTransport(mcParticle& p) override;
	void beginTransportInside(mcParticle& p) override;

	double getDistanceOutside(mcParticle& p) const override;

	mcTransport* getInternalTransportByName(const char* name) override;

	// ���������� ������� ���������� �������.
	// ����������� ���������� ���������� ���������� ��� ����������� 
	// � ����� ��������� ������� �������, ����� ��� ��������� � ������� �������.
	void addTransport(mcTransport* t);

	// ���������� �� ��������� �������� ��� ����������� ��������� ���� ������
	void completeInit();

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	std::vector<mcTransport*> chainTransports_;
};
