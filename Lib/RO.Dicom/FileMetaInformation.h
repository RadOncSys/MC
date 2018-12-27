// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"
#include "DataSet.h"

namespace RO
{
	namespace Dicom
	{
		class RODICOM_API FileMetaInformation
		{
		public:
			FileMetaInformation();

			DataSet& Elements() { return elements_; }

			void ReadFromStream(std::istream& is);
			void WriteToStream(std::ostream& os);

			static std::shared_ptr<FileMetaInformation> CreateDefault();
			std::shared_ptr<TransferSyntaxItem> Ts();

		private:
			DataSet elements_;
			std::shared_ptr<TransferSyntaxItem> ts_;
		};
	}
}
