// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include <memory>
#include <string>

namespace RO
{
	namespace Dicom
	{
		class RODICOM_API TransferSyntaxItem
		{
		public:
			TransferSyntaxItem(const char* name, const char* uid, bool bLittleEndian, bool bEncapsulated, bool bExplicitVr, bool bDeflate, bool bLossy, bool bLossless) :
				LittleEndian(bLittleEndian), Encapsulated(bEncapsulated), ExplicitVr(bExplicitVr), Deflate(bDeflate), Lossy(bLossy),
				Lossless(bLossless), Name(name), Uid(uid) {}

			const bool LittleEndian;
			const bool Encapsulated;
			const bool ExplicitVr;
			const bool Deflate;
			const bool Lossy;
			const bool Lossless;
			const std::string Name;
			const std::string Uid;
		};

		namespace TransferSyntax
		{
			RODICOM_API std::shared_ptr<TransferSyntaxItem> Item(const char* uid);
			RODICOM_API void Populate();
		}

		typedef std::shared_ptr<TransferSyntaxItem> TransferSyntaxType;
	}
}
