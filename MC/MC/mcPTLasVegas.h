// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"

// ����� ���������� � ������� Las Vegas, ������������ ��� ������������ ���������� ������.
// � ������� ��������� ������� ����� �������������, � ������� ����������� ���������, 
// ����� � ��������� (��������� Z=0).
class mcPTLasVegas : public mcTransport
{
public:
	mcPTLasVegas(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x);
	virtual ~mcPTLasVegas(void);

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;
	double getDNearInside(const geomVector3D& p) const override;

	// ����� ��� ������������
	double A() const { return a_; }
	double Z() const { return z_; }
	const std::vector<double>& Ds() const { return ds_; }
	const std::vector<double>& Hs() const { return hs_; }
	const std::vector<double>& Xs() const { return xs_; }
	const std::vector<double>& Ys() const { return ys_; }

protected:
	// ������ ����������� ��� ��������� ����� �������������� ��������� (5�5)
	double a_;	// ������� ����������� �������
	double z_;	// ������� �������

	// �������� ���������
	std::vector<double> ds_;

	// ������� ���������
	std::vector<double> hs_;

	// ������� ����� ����� ���������
	std::vector<double> xs_;
	std::vector<double> ys_;

	// ��������� ����������, ������������ ��� ����� ���������� ������������ ���������
	double x1_, x2_, y1_, y2_, z2_;
	double ax_, ay_;	// ������������ �������������, ������������ ������� ���������
};
