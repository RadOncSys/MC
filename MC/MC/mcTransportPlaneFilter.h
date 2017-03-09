// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// ����� ������������ �������.
// ���� - ����������� ������ ����� ����������� ��������� 
// (�������� ���������������� ��� Z) � ����� ���������� � ���.
// � ����������� �� ������������� ������ �� ���������� � ����������� ���������� 
// �������� ������������ � �������� ����������� ������.
class mcTransportPlaneFilter : public mcTransport
{
public:
	mcTransportPlaneFilter(void);
	mcTransportPlaneFilter(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcTransportPlaneFilter(void);

	// ������ ���������� ����������, ��� ��� �� ����� ��������� ��������� ������.
	// �� ������ ����� ���������������� � ���������� ������.
	void beginTransport(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;
};
