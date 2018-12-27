// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include <string>

namespace RO
{
	namespace Dicom
	{
		//const std::string ImplementationClassUID = "1.2.840.416.480.5705.1234";//pretty much random, no policy implemented.
		const std::string UidImplementationClassUID = "1.2.826.0.1.3680043.2.1553";//UID prefix supplied by Medical Connection in Sep 2006
		const std::string UidImplementationVersionName = "DICOMLIB2008";
	}
}
