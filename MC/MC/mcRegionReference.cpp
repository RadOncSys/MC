#include "mcRegionReference.h"

mcRegionReference::mcRegionReference(void) : medidx_(0), subidx_(0), idx_(0)
{
	gidx_[0] = 0;
	gidx_[1] = 0;
	gidx_[2] = 0;
}

mcRegionReference::~mcRegionReference(void)
{
}
