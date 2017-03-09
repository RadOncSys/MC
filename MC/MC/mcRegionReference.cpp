#include "mcRegionReference.h"

mcRegionReference::mcRegionReference(void)
	:type_(MCRT_LIN_IDX)
	, medidx_(0)
	, idx_(0)
{
	gidx_[0] = 0;
	gidx_[1] = 0;
	gidx_[2] = 0;
}

mcRegionReference::~mcRegionReference(void)
{
}
