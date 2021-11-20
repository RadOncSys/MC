// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2021] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include "mcMedium.h"
#include "mcDefs.h"

// ����� �������� ���������� ���������� ����� ��� ���������� ��������.
class mcMediumNeutron : public mcMedium
{
public:
	mcMediumNeutron(void);
	virtual ~mcMediumNeutron(void);

	virtual void read(istream& is);
};
