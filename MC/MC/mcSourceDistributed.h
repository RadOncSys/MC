// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"

// ���� ������������� ����������
enum mc_distr_t { MCP_CONST_DENS = 0, MCP_GAUSSIAN, MCP_WATERBAG, MCP_K_V, MCP_ARBITRARY };


class mcSourceDistributed : public mcSource
{
public:
	mcSourceDistributed(void);
	mcSourceDistributed(const char* name, int nThreads,
		mc_particle_t type, double ke, const geomVector3D& p, const geomVector3D& v,
		mc_distr_t distr, double Emitx, double rbeamx, double beamAnglex, double Emity,
		double rbeamy, double beamAngley);

	void init(mc_particle_t type       // ��� ������
		, double ke               // ������������ �������
		, const geomVector3D& p   // ����� �������� ������
		, const geomVector3D& v   // ����������� �������� (��������� ������)
		, mc_distr_t distr        // ��� �������������
		, double Emitx  // �-��������, � ��. "���������������", ��� ��, �.�. ������ rbeam*rbeta, rbeta = vrmax/C, C - �������� �����
		, double rbeamx // ������� �� �
		, double beamAnglex // ������� ����������� �� � (��-�� ������, �� ������� � ����������)
		, double Emity            //
		, double rbeamy           //  ��� �� �� ����� ��� y
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
	mc_distr_t distr_;  // ��� �������������
	double Emitx_;  // �-��������, � ��. "���������������", ��� ��, �.�. ������ rbeam*rbeta, rbeta = vrmax/C, C - �������� �����
	double rbeamx_; // ������� �� �
	double rbetax_; // Tmitx/rbeamx
	double beamAnglex_; // ������� ����������� �� � (��-�� ������, �� ������� � ����������)
	double Emity_;            //
	double rbeamy_;           //  ��� �� �� ����� ��� y
	double rbetay_;           //
	double beamAngley_;       //
	double beta0_;   // Vz/C, sqrt(ke/Ep/(ke/Ep+1)) - � ������ �����������

	mc_particle_t type_;
	double ke_;       // � ���                (���������!!!!!!!!!!!!!!!!!!!)
	geomVector3D p_;  // ���������� ������ ���������
	geomVector3D v_;  // ���������� ��������, ���� ��������� ������ (0,0,1)
	int q_; // �����
	// ��� ������������� ������������� - ��������� �� ������� �������, ������������ double �� 0 �� 1
	//  double (*fdistr_)(double x,double y,double px,double py);
};
