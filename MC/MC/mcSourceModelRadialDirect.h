// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcsource.h"
#include <vector>

/// <summary>
/// ����� ������ ��������� ��������� � ���������� ����������, ����������� ������� ����� ������ (�������� ������������).
/// ������ � ������������� ������ (r, Ux, Uy, E), ��� Ux � Uy - �������� ����������� �������� � ����������
/// ����������� ����� ��� ��������� � ������� ��������� ����� � ���������������� ��.
/// �����, ������������� � ����� ������ ������� ��� ����� � ���������������� ����������.
/// ��� ����� �������� �� ������ ������ ����������� 8 �� ������.
/// ��������, ���������� ��������� � ������������� ����������� � ��������� ������ ��� ����� ����� ���-�� ���� 
/// 100 ��������� ������ � ���������� ����� ������������
/// </summary>
class mcSourceModelRadialDirect : public mcSource
{
public:
	// name		- ������ �������� �������
	// nThreads - ���������� ������������� ������� (��� ������ ���� ������ � ����� ������� ������).
	// z0		- ��������� ��������� ������������ �������, � �������� �� ��������
	mcSourceModelRadialDirect(const char* name, int nThreads, double z0);
	~mcSourceModelRadialDirect();

	void sample(mcParticle& p, mcThread* thread) override;

	// ������������ ������ (��� �� ���������)
	void Init(double maxEnergy, double maxR);

	// �������� ��������������� ����� ������ ��������� ������.
	// ��������� ��������� ������� �������� � ������.
	void SetSplitting(unsigned nSplits);

	// ���������� ���������� ������ �� ������� ���������.
	void SetParticles(unsigned short* data, unsigned nThreadSize, const std::vector<unsigned>& indexes);
	
	// ������ � �������� ��� ��������� �����
	const unsigned short* GetParticles() const { return particles_; }
	unsigned GetNoofParticles() const { return noofParticles_; }

	/// <summary>
	/// �������� ������ ������ � ������. ������ ������������� ������ �������.
	/// ������� ���������� ������� ������� ���������� �� � ������� free();
	/// �������� size �������� ����� ���������� ������ � ������
	/// </summary>
	void* saveToMemory(int& size);

	/// <summary>
	/// �������������� ������ �� ������ ������.
	/// </summary>
	void readFromMemory(void* buffer);

	void dumpVRML(ostream& os) const override;

	friend std::ostream& operator << (std::ostream&, const mcSourceModelRadialDirect&);

protected:
	// ��������� ��������� ��������� �� ��� �������
	double z0_;

	// ������������, �� ������� ����� ��������, ����� �������� ���������� �������� ����������
	double energyScale_;
	double rScale_;				// �������� �������������� ��������
	double radialAngleScale_;	// ���������� Ux �������� �������� ([-1,1] -> [0,2])
	double azimutAngleScale_;	// ���������� Uy �������� �������� ([-1,1] -> [0,2])

	// ���������� ������ � ����� ����������
	unsigned noofParticles_;

	// ���������� ������ ��� �����������
	unsigned noffsplits_;

	// ������ ������� ������� ��� ���������
	std::vector<unsigned> currentParticleIdx_;

	// ������ ������� ������� �� ���������, ������� ������������� �������������� ����������.
	// ���������� ���� ������ ������ �� ��� ������, �� � ������� ������ ���� �������.
	std::vector<unsigned> currentSplitIdx_;

	// �������� ���� � �������� ��� ���������� ������
	std::vector<double> splitStartAngle_;
	
	// ��������� �� ������� ����������� ��� �� ����
	double da_;

	// ����������� ���� ������������ ������
	double dw_;

	// ������ ���� ������ � ��� ����� � � ������ ��������� ���������.
	// � ��������� ������ �� ��������� ����� �� ������ ����� ��� ������� ������.
	// ����� ��������� �������� � �������� ������ ������
	unsigned short* particles_;
};
