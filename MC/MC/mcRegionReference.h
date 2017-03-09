// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

enum mc_region_idx_t { MCRT_LIN_IDX = 0, MCRT_3D_IDX, MCRT_BOTH_IDX };

// Класс, отвечающий за обозначение региона.
// Это может быть линейный индекс региона 
// или 3 индекса 3D сетки, или и то и другое
class mcRegionReference
{
public:
	mcRegionReference(void);
	~mcRegionReference(void);

	enum mc_region_idx_t type_;
	short medidx_;  // индекс среды
	int idx_;      // индекс региона
	short gidx_[3];
};
