// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

// Класс, отвечающий за обозначение региона.
// Это может быть линейный индекс региона 
// или 3 индекса 3D сетки, или и то и другое
class mcRegionReference
{
public:
	mcRegionReference(void);
	~mcRegionReference(void);

	short medidx_;  // индекс среды
	short subidx_;  // дополнительный индекс, например, внутри тела в ячейке или снруже в случае гребенчатого фильтра
	int idx_;       // индекс региона
	short gidx_[3];
};
