// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

enum mc_region_idx_t { MCRT_LIN_IDX = 0, MCRT_3D_IDX, MCRT_BOTH_IDX };

// �����, ���������� �� ����������� �������.
// ��� ����� ���� �������� ������ ������� 
// ��� 3 ������� 3D �����, ��� � �� � ������
class mcRegionReference
{
public:
	mcRegionReference(void);
	~mcRegionReference(void);

	enum mc_region_idx_t type_;
	short medidx_;  // ������ �����
	int idx_;      // ������ �������
	short gidx_[3];
};
