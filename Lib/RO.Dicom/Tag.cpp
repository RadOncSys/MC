#include "pch.h"
#include "Tag.h"

using namespace RO::Dicom;

UINT32 RO::Dicom::MakeTag(UINT16 group, UINT16 item) { return (UINT32(group) << 16) | item; };
UINT16 RO::Dicom::TagGroup(UINT32 tag) { return UINT16(tag >> 16); };
UINT16 RO::Dicom::TagElement(UINT32 tag) { return UINT16(tag); };
