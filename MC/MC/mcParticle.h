// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcRegionReference.h"
#include "../geometry/vec3d.h"

class mcScoreTrack;
class mcThread;
class mcTransport;

// ���� ������ ������� � �������� �������� �������� ����������
enum mc_particle_t { MCP_PHOTON = 0, MCP_NEGATRON, MCP_POSITRON, MCP_PROTON, MCP_NEUTRON, MCP_NTYPES };

class mcParticle
{
public:
	mcParticle(void);
	mcParticle(mc_particle_t pt, int pq, double pke, const geomVector3D& pp, const geomVector3D& pu);
	mcParticle(const mcParticle& p);
	~mcParticle(void);

public:
	enum mc_particle_t t;
	int q;
	double ke;
	geomVector3D p;
	geomVector3D plast;	// ����� ���������� ��������������
	geomVector3D u;
	double dnear;
	mcRegionReference region;
	double weight;

	// ���������� �� ���������� �������������� � �������� ���� �������� �������
	double mfps;

	// ������������� ��������� �������, � ������� ��������� ������� � ������ ������
	double regDensityRatio;

	// ��������� �� ������������ ������,
	// � ������� ������� ��������� � ������ ������
	mcTransport* transport_;

	// ��������� �� ������������ ������ � ����������� ����������� � ������ ����������� 
	// �������� ���������� ������ ������
	mcTransport* transportNearest_;

	// ����� �������� �������
	//mcRegionReference reg_born;

	// ����� ���������� ��������������
	//mcRegionReference reg_last_scatter;

	// ��������� �� ������, ���������� ����������� ������ ����� ��������
	mcScoreTrack* trackScore_;

	// ��������� �� ������, ���������� ����������� ��� ������ ������.
	// ��������� ������ ��� �������� ������������ �����������,
	// ����������� � ������������� ����������� ��� ������ ������,
	// �� ������ ����� ��� ����������� ������� ������ ���������� ��� ����������, 
	// ������� ���������� ��������� ����� � ������������� ������.
	mcThread* thread_;

	// ������� ���������� ����� ������� ������ � ������ 
	// �������� �������. ���������� ������� �������� �� ��������
	int regionBirth;

	// ������� ���������� ����� ������� ������ � ������ 
	// BeginTransport ���������, ��� ������� � ��� ����������.
	int regionFlags;

	// ��� ������������ ����������� (������� ��� ���������) � embedded transport.
	// �� ������������ �������� ������� ��� �������� �� ���, ��� ��������������� 
	// ������� ���������� ���������� ��������� ������������� ��������� ������� �����.
	enum temb_shit_t : short { Undefined = 0, External, Internal };
	temb_shit_t exitSurface_;
};
