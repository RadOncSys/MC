// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once

#ifdef RODICOM_EXPORTS
#define RODICOM_API __declspec(dllexport)
#else
#define RODICOM_API __declspec(dllimport)
#endif

#include <Windows.h>

#pragma warning(disable:4251)  // Disable warning about needs to have dll-interface to be used by clients
