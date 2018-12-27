// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include <iostream>
#include <memory>

namespace RO
{
	namespace Dicom
	{
		class DataSet;
		class FileMetaInformation;

		class RODICOM_API DicomObject
		{
		public:
			DicomObject(void);
			~DicomObject(void);

			std::shared_ptr<DataSet> Dset;
			std::shared_ptr<FileMetaInformation> MetaInfo;

			void ReadFromStream(std::istream& is, bool skipPixels = false);
			void ReadFromFile(std::string fname, bool skipPixels = false);
			void ReadFromBuffer(const char* data, unsigned int len);

			// Запись данных включая мета-данные
			void WtiteToStream(std::ostream& os);
			void WtiteToFile(std::string fname);

			void* PixelData();
			long PixelDataLength();	// in bytes
		};
	}
}
