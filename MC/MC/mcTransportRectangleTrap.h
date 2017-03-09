// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransport.h"

// ����� ���������� ������� ��������� ��� �������, �� ����������� ���
// ������� ����� � ������������� ����������� � �������� � �������� �������������.
class mcTransportRectangleTrap : public mcTransport
{
public:
	mcTransportRectangleTrap(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, double width, double length);
	~mcTransportRectangleTrap(void);

	// ������ ���������� ����������, ��� ��� �� ����� ��������� ��������� ������.
	// �� ������ ����� �������������.
	void beginTransport(mcParticle& p) override;

	void dumpVRML(ostream& os)const override;

protected:
	double swidth_;
	double slength_;
};
