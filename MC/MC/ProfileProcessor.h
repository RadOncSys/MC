// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <vector>
#include <memory>

// ����� ������� ������� ��������.
class ProfileProcessor
{
public:
	/// <summary>
	/// ����������� ���������� ����������� �������.
	/// </summary>
	/// <param name="x">���������� ����� �������</param>
	/// <param name="d">�������� ���</param>
	/// <param name="width">����������� ������ ������� �� ����������</param>
	/// <param name="penumbra">����������� ������ ���������</param>
	/// <param name="penumbra">��������������� ������, ��������������� ��������� ������� �� 80% ���� 3 ��</param>
	static void ProfileParameters(const std::vector<double>& x, const std::vector<double>& d,
		std::vector<double>& dd, double& width, double& penumbra, double& dinde);

	/// <summary>
	/// ����������� ������� ������������� ������� ����� � FanLine RZ ���������.
	/// ���� ����������� � �� ������� � ����������� ����������� ��� ������������� ���������.
	/// ����������� ����� �� ���������� �� ����� �������� ������� ���� � �� ����� � ��� ������� ��� 1.5 ��.
	/// ��� ���� � ���� ����� ������������ �������������� ����������� ����������� ���������������� ������ ���������.
	/// �� ������� ������������� �� ����������� ������� build-up �������� ��������.
	/// ������� ���������� �� ���������� Savitzky�Golay ������� �������.
	/// </summary>
	/// <param name="srcMatrix">�������� ������ ������� ����� (������ ��������� �� ������� ������� �� ��������)</param>
	/// <param name="nr">������ ������� �� �������</param>
	/// <param name="nz">������ ������� �� �������</param>
	/// <param name="rstep">��� �� ������� � �� (��� ����������� �������� ����������� �� �������)</param>
	static std::shared_ptr<std::vector<std::vector<double>>>
		SmoothFanRZ(const std::vector<std::vector<double>>& srcMatrix, unsigned nr, unsigned nz, double rstep);

	// ����������� ���������� ������� ��������� ������������� ������� ���������
	// nr - ����������� ����� �� ������, ������������ �����������.
	static void SmoothRZProfile(std::vector<double>&, unsigned nr);

	// ���������� ����������� Savitzky�Golay � ���������� ���������.
	static void SmoothSG1D(std::vector<double>& p, unsigned m, unsigned nl, unsigned nr);

	// ��������� ����������� Savitzky�Golay.
	static std::unique_ptr<std::vector<std::vector<double>>> SmoothSG2D(const std::vector<std::vector<double>>& data);
};
