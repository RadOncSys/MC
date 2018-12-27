// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once
#include "RO.Dicom.h"

#include "DicomObject.h"
#include <vector>

namespace RO
{
	namespace Dicom
	{
		class dcmImage;
		struct ddr_patient_record;

		struct RODICOM_API ddr_seriesleaf_record
		{
			ddr_seriesleaf_record(){}
			std::string uid_;
			std::string file_;
			std::shared_ptr<DicomObject> dcmobj;
		};

		struct RODICOM_API ddr_series_record
		{
			std::string uid_;
			std::string number_;
			std::string datetime_;
			std::string modality_;
			std::string descr_;
			std::shared_ptr<std::vector<ddr_seriesleaf_record>> leafs_;
		};

		struct RODICOM_API ddr_study_record
		{
			ddr_study_record(){}
			std::string uid_;
			std::string id_;
			std::string institution_;
			std::string datetime_;
			std::string number_;
			std::string descr_;
			std::shared_ptr<ddr_patient_record> patient_;
			std::shared_ptr<std::vector<ddr_series_record>> series_;
		};

		struct RODICOM_API ddr_patient_record
		{
			std::string id_;
			std::string name_;
			std::string birthDate_;
			std::string sex_;
			std::shared_ptr<std::vector<ddr_study_record>> studies_;
		};

		///  ласс обслуживани€ каталога исследований включа€ списки изображений
		class RODICOM_API DicomDir
		{
		public:
			DicomDir(std::string dcmdirfile);
			~DicomDir(void);

			void GetPatientList(std::shared_ptr<std::vector<ddr_patient_record>>);
			void GetStudyList(std::shared_ptr<std::vector<ddr_study_record>>, const std::string& patientID);
			void GetSeriesList(std::shared_ptr<std::vector<ddr_series_record>>, 
				const std::string& patientID, const std::string& studyInstanceUID);
			void GetSeriesLeafList(std::shared_ptr<std::vector<ddr_seriesleaf_record>>,
				const std::string& patientID, const std::string& studyInstanceUID, const std::string& seriesInstanceUID);
			void CreateCatalog(std::shared_ptr<std::vector<ddr_patient_record>>, const std::string& modality);

			/// —канирование файлов в дереве директорий и генераци€ каталога
			static void CreateCatalogFromDcmFileTree(const std::string& rootdir, 
				std::shared_ptr<std::vector<std::string>> flist, const std::string& modality,
				std::vector<ddr_patient_record>& catalog, const std::string& studyUid, const std::string& seriesUid);

			/// —канирование списка изображений и генераци€ каталога
			static void createCatalogFromDcmObjList(std::shared_ptr<std::vector<std::shared_ptr<DicomObject>>> olist,
				const std::string& modality, std::shared_ptr<std::vector<ddr_patient_record>> catalog);

		protected:
			/// ƒобавление в каталог новой записи об изображении
			static std::shared_ptr<ddr_seriesleaf_record> addRecordFromImage(
				std::shared_ptr<std::vector<ddr_patient_record>> catalog, std::shared_ptr<DicomObject> img);

			/// —оздание каталога только с теми данными, которые относ€тс€ к указанной модальности
			static void copyCatalogWithModalityFilter(std::shared_ptr<std::vector<ddr_patient_record>> toCatalog,
				std::shared_ptr<std::vector<ddr_patient_record>> fromCatalog, const std::string& modality);

		protected:
			std::string ddrfname_;
			void *object_; // DICOMDIR file image
		};
	}
}
