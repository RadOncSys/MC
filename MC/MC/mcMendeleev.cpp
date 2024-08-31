#include "mcMendeleev.h"

mcMendeleev::mcMendeleev()
{
	IsNecessary.resize(119, false);
	IsLoad.resize(119, false);
}

mcMendeleev::mcMendeleev(const mcMendeleev& t)
{
	IsNecessary.resize(119, false);
	IsLoad.resize(119, false);
	for (int i = 0; i < IsNecessary.size(); i++)
	{
		IsNecessary[i] = t.IsNecessary[i];
		IsLoad[i] = t.IsLoad[i];
	}
}
