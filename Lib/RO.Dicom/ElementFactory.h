// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include "Element.h"
#include <memory>

namespace RO
{
	namespace Dicom
	{
		class RODICOM_API ElementFactory
		{
		public:
		    static std::shared_ptr<Element> MakeElement(UINT32 tag, Vr::VR vr);
		};
	}
}
