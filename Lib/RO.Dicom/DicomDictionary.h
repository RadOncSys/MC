// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include "Vr.h"
#include <memory>
#include <string>

namespace RO
{
	namespace Dicom
	{
		class DicomDictionaryItem
		{
		public:
			DicomDictionaryItem(UINT32 tag, enum Vr::VR vr, std::string vm, std::string keyword, std::string commnent) :
				Tag(tag), Vr(vr), Vm(vm), Keyword(keyword), Commnent(commnent) {}
			DicomDictionaryItem(UINT32 tag, enum Vr::VR vr, std::string vm) :
				Tag(tag), Vr(vr), Vm(vm) {}

			const UINT32 Tag;
			const enum Vr::VR Vr;
			const std::string Vm;
			const std::string Keyword;
			const std::string Commnent;
		};

		namespace DicomDictionary
		{
			RODICOM_API std::shared_ptr<DicomDictionaryItem> Item(UINT32 tag);
			RODICOM_API std::shared_ptr<DicomDictionaryItem> Item(const std::string& keyword);
			RODICOM_API std::shared_ptr<DicomDictionaryItem> ItemByKey(const char* keyword);
			RODICOM_API void Populate();
			RODICOM_API void Populate(const char* data, unsigned int len);
			RODICOM_API void Populate(std::istream& is);

			RODICOM_API Vr::VR StringToVr(const std::string& s);
			RODICOM_API std::string VrToString(Vr::VR vr);
		}
	}
}
