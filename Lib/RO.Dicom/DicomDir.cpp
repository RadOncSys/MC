#include "pch.h"
#include "DicomDir.h"
//
//
//DicomDir::DicomDir(void)
//{
//}
//
//
//DicomDir::~DicomDir(void)
//{
//}
//
//
//
//// ddr.cpp:  Вспомогательные функции для инкапсуляции обращений к библиотеке CTN
//// Проект:	 AMPhora
//// Автор:		 Горлачев Г.Е. (ggorl@dol.ru)
///////////////////////////////////////////////////////////////////////////////
//
//#include "ddr.h"
//#include "dcmImage.h"
//#include "dcmelement.h"
//#include "dcmaux.h"
//#include "../util/datetime.h"
//#include "../util/file.h"
//
//dcmuDICOMDIR::dcmuDICOMDIR(const char* fname)
//:ddrfname_(fname)
//{
//	//CONDITION cond = DCM_OpenFile(fname, DCM_ORDERLITTLEENDIAN, (DCM_OBJECT**)&object_);
//	//Could not open dicom file as expected.  Trying Part 10 format.
//	//if(cond != DCM_NORMAL)
//	CONDITION cond = DCM_OpenFile(fname, DCM_PART10FILE, (DCM_OBJECT**)&object_);
//	if (cond != DCM_NORMAL)
//    dcmuRiseException(cond);
//}
//
//dcmuDICOMDIR::~dcmuDICOMDIR()
//{
//	if(object_ != NULL)
//		DCM_CloseObject((DCM_OBJECT**)&object_);
//}
//
//void dcmuDICOMDIR::getPatientList(vector<ddr_patient_record>& rset)
//{
//  LST_HEAD* patientList = LST_Create();
//  CONDITION cond = DDR_GetPatientList((DCM_OBJECT**)&object_, &patientList);
//  unsigned i,n=LST_Count(&patientList);
//  rset.resize(n);
//  DDR_PATIENT* patientNode = (DDR_PATIENT*)LST_Dequeue(&patientList);
//  for(i=0; i<n && patientNode != NULL; i++){
//    ddr_patient_record& r=rset[i];
//    r.id_=patientNode->PatientID;
//    r.name_=patientNode->PatientName;
//    r.birthDate_=patientNode->BirthDate;
//    r.sex_=patientNode->Sex;
//    free(patientNode);
//    patientNode = (DDR_PATIENT*)LST_Dequeue(&patientList);
//  }
//  LST_Destroy(&patientList);
//}
//
//void dcmuDICOMDIR::getStudyList(vector<ddr_study_record>& rset
//                                ,const char* patientID)
//{
//  LST_HEAD* studyList = LST_Create();
//  CONDITION cond = DDR_GetStudyList((DCM_OBJECT**)&object_, patientID, &studyList);
//  unsigned i,n=LST_Count(&studyList);
//  rset.resize(n);
//  DDR_STUDY* studyNode = (DDR_STUDY*)LST_Dequeue(&studyList);
//  for(i=0; i<n && studyNode != NULL; i++){
//    ddr_study_record& r=rset[i];
//    r.uid_=studyNode->StudyInstanceUID;
//    r.id_=studyNode->StudyID;
//    rtuW32DateConverter dt;
//    dt.setDateFromDicomString(studyNode->StudyDate);
//    dt.setTimeFromDicomString(studyNode->StudyTime);
//    r.datetime_=dt.getDateTimeString();
//    r.number_=studyNode->AccessionNumber;
//    r.descr_=studyNode->StudyDescription;
//    free(studyNode);
//    studyNode = (DDR_STUDY*)LST_Dequeue(&studyList);
//  }
//  LST_Destroy(&studyList);
//}
//
//void dcmuDICOMDIR::getSeriesList(vector<ddr_series_record>& rset
//                                 ,const char* patientID
//                                 ,const char* studyInstanceUID)
//{
//  LST_HEAD* seriesList = LST_Create();
//  CONDITION cond = DDR_GetSeriesList((DCM_OBJECT**)&object_, patientID, studyInstanceUID, &seriesList);
//  unsigned i,n=LST_Count(&seriesList);
//  rset.resize(n);
//  DDR_SERIES* seriesNode = (DDR_SERIES*)LST_Dequeue(&seriesList);
//  for(i=0; i<n && seriesNode != NULL; i++){
//    ddr_series_record& r=rset[i];
//    r.uid_=seriesNode->SeriesInstanceUID;
//    r.number_=seriesNode->SeriesNumber;
//    r.modality_=seriesNode->Modality;
//    free(seriesNode);
//    seriesNode = (DDR_SERIES*)LST_Dequeue(&seriesList);
//  }
//  LST_Destroy(&seriesList);
//}
//
//void dcmuDICOMDIR::getSeriesLeafList(vector<ddr_seriesleaf_record>& rset
//                                     ,const char* patientID
//                                     ,const char* studyInstanceUID
//                                     ,const char* seriesInstanceUID)
//{
//  LST_HEAD* leafList = LST_Create();
//  CONDITION cond = DDR_GetSeriesLeafList((DCM_OBJECT**)&object_, patientID, studyInstanceUID, seriesInstanceUID, &leafList);
//  unsigned i,n=LST_Count(&leafList);
//  rset.resize(n);
//  DDR_SERIES_LEAF* leafNode = (DDR_SERIES_LEAF*)LST_Dequeue(&leafList);
//  for(i=0; i<n && leafNode != NULL; i++){
//    ddr_seriesleaf_record& r=rset[i];
//    r.file_=leafNode->FileID;
//    r.uid_=leafNode->SOPInstanceUID;
//    free(leafNode);
//    leafNode = (DDR_SERIES_LEAF*)LST_Dequeue(&leafList);
//  }
//  LST_Destroy(&leafList);
//}
//
//void dcmuDICOMDIR::createCatalog(vector<ddr_patient_record>& catalog, const char* modality)
//{
//  unsigned i,j,k;
//  vector<ddr_patient_record> catalogTmp;
//  getPatientList(catalogTmp);
//  for(i=0;i<catalogTmp.size();i++){
//    ddr_patient_record& patient = catalogTmp[i];
//    getStudyList(patient.studies_
//                 ,patient.id_.c_str());
//
//    for(j=0;j<patient.studies_.size();j++){
//      ddr_study_record& study = patient.studies_[j];
//      study.patient_ = &patient;
//      getSeriesList(study.series_
//                    ,patient.id_.c_str()
//                    ,study.uid_.c_str());
//
//      for(k=0;k<study.series_.size();k++){
//        ddr_series_record& series = study.series_[k];
//        getSeriesLeafList(series.leafs_
//                          ,patient.id_.c_str()
//                          ,study.uid_.c_str()
//                          ,series.uid_.c_str());
//        try {
//          if(!series.leafs_.empty()) {
//            string rootdir = FilePath(ddrfname_);
//            dcmImage img;
//            img.read((rootdir+"\\"+series.leafs_[0].file_).c_str());
//
//            study.institution_ = img.elementString(DCM_IDINSTITUTIONNAME);
//            study.institution_ += " / ";
//            study.institution_ = img.elementString(DCM_IDMANUFACTURERMODEL);
//            study.institution_ += " / ";
//            study.institution_ = img.elementString(DCM_IDSTATIONNAME);
//
//            rtuW32DateConverter dt;
//            dt.setDateFromDicomString(img.elementString(DCM_IDSERIESDATE));
//            dt.setTimeFromDicomString(img.elementString(DCM_IDSERIESTIME));
//            series.datetime_=dt.getDateTimeString();
//
//            series.descr_=img.elementString(DCM_IDSERIESDESCR);
//          }
//        }
//        catch(...){}
//      }
//    }
//  }
//  copyCatalogWithModalityFilter(catalog, catalogTmp, modality);
//}
//
////
//// добавление в каталог новой записи об изображении
////
//ddr_seriesleaf_record& 
//dcmuDICOMDIR::addRecordFromImage(vector<ddr_patient_record>& catalog, dcmImage& img)
//{
//  unsigned j;
//  ddr_patient_record precord;
//  precord.id_=img.elementString(DCM_PATID);
//  precord.name_=img.elementString(DCM_PATNAME);
//  precord.birthDate_=img.elementString(DCM_PATBIRTHDATE);
//  precord.sex_=img.elementString(DCM_PATSEX);
//  // поиск
//  ddr_patient_record* patient=0;
//  for(j=0;j<catalog.size();j++){
//    ddr_patient_record& p=catalog[j];
//    if(p.id_==precord.id_ && p.name_==precord.name_){
//      patient=&p;
//      break;
//    }
//  }
//  if(!patient){
//    catalog.push_back(precord);
//    patient=&catalog.back();
//  }
//
//  // исследование
//  ddr_study_record srecord;
//  srecord.uid_=img.elementString(DCM_RELSTUDYINSTANCEUID);
//  srecord.id_=img.elementString(DCM_RELSTUDYID);
//
//  srecord.institution_ = img.elementString(DCM_IDINSTITUTIONNAME);
//  srecord.institution_ += " / ";
//  srecord.institution_ = img.elementString(DCM_IDMANUFACTURERMODEL);
//  srecord.institution_ += " / ";
//  srecord.institution_ = img.elementString(DCM_IDSTATIONNAME);
//
//  rtuW32DateConverter dt;
//  dt.setDateFromDicomString(img.elementString(DCM_IDSTUDYDATE));
//  dt.setTimeFromDicomString(img.elementString(DCM_IDSTUDYTIME));
//  srecord.datetime_=dt.getDateTimeString();
//
//  srecord.number_=img.elementString(DCM_IDACCESSIONNUMBER);
//  srecord.descr_=img.elementString(DCM_IDSTUDYDESCRIPTION);
//  // поиск
//  ddr_study_record* study=0;
//  for(j=0;j<patient->studies_.size();j++){
//    ddr_study_record& s=patient->studies_[j];
//    if(s.uid_==srecord.uid_){
//      study=&s;
//      break;
//    }
//  }
//  if(!study){
//    patient->studies_.push_back(srecord);
//    study=&patient->studies_.back();
//    study->patient_ = patient;
//  }
//
//  // серия
//  ddr_series_record serrecord;
//  serrecord.uid_=img.elementString(DCM_RELSERIESINSTANCEUID);
//  serrecord.number_=img.elementString(DCM_RELSERIESNUMBER);
//  serrecord.modality_=img.elementString(DCM_IDMODALITY);
//
//  dt.setDateFromDicomString(img.elementString(DCM_IDSERIESDATE));
//  dt.setTimeFromDicomString(img.elementString(DCM_IDSERIESTIME));
//  serrecord.datetime_=dt.getDateTimeString();
//
//  serrecord.descr_=img.elementString(DCM_IDSERIESDESCR);
//  // поиск
//  ddr_series_record* series=0;
//  for(j=0;j<study->series_.size();j++){
//    ddr_series_record& s=study->series_[j];
//    if(s.uid_==serrecord.uid_){
//      series=&s;
//      break;
//    }
//  }
//  if(!series){
//    study->series_.push_back(serrecord);
//    series=&study->series_.back();
//  }
//
//  // файл
//  ddr_seriesleaf_record leaf;
//  series->leafs_.push_back(leaf);
//
//  return series->leafs_.back();
//}
//
////
//// сканирование списка изображений и генерация каталога
////
//void dcmuDICOMDIR::createCatalogFromDcmFileTree(const string& rootdir
//                                                ,const vector<string>& flist
//                                                ,const char* modality
//                                                ,vector<ddr_patient_record>& catalog
//                                                ,const char* studyUid
//                                                ,const char* seriesUid)
//{
//  string modlist(modality);
//  for(unsigned i=0;i<flist.size();i++){
//    try{
//      dcmImage img;
//      img.read((rootdir+"\\"+flist[i]).c_str());
//      if(modality && modlist.find(img.elementString(DCM_IDMODALITY))==string::npos) 
//        continue;
//      if(studyUid && strcmp(studyUid, img.elementString(DCM_RELSTUDYINSTANCEUID)) != 0)
//        continue;
//      if(seriesUid && strcmp(seriesUid, img.elementString(DCM_RELSERIESINSTANCEUID)) != 0)
//        continue;
//      ddr_seriesleaf_record& leaf = addRecordFromImage(catalog, img);
//      leaf.uid_=img.elementString(DCM_IDSOPINSTANCEUID);
//      leaf.file_=flist[i];
//    }
//    catch(...){}
//  }
//}
//
////
//// сканирование списка изображений и генерация каталога
////
//void dcmuDICOMDIR::createCatalogFromDcmObjList(const vector<dcmImage*>& olist
//                                               ,const char* modality
//                                               ,vector<ddr_patient_record>& catalog)
//{
//  for(unsigned i=0;i<olist.size();i++){
//    dcmImage& img = *olist[i];
//    if(modality && strcmp(modality, img.elementString(DCM_IDMODALITY)) !=0 ) 
//      continue;
//    ddr_seriesleaf_record& leaf = addRecordFromImage(catalog, img);
//    leaf.uid_=img.elementString(DCM_IDSOPINSTANCEUID);
//    // файл
//    leaf.ctx_=(void*)&img;
//  }
//}
//
//// создание каталога только с теми данными, которые относятся к указанной модальности
//void dcmuDICOMDIR::copyCatalogWithModalityFilter(vector<ddr_patient_record>& toCatalog
//                                                 ,const vector<ddr_patient_record>& fromCatalog
//                                                 ,const char* modality)
//{
//  if(modality == 0) {
//    toCatalog = fromCatalog;
//    return;
//  }
//  string modlist(modality);
//
//  unsigned i,j,k;
//  for(i=0;i<fromCatalog.size();i++)
//  {
//    const ddr_patient_record& patient = fromCatalog[i];
//
//    // Проверка наличия нужных серий у пациента
//    bool isPatientValid = false;
//    for(j=0; j < patient.studies_.size(); j++)
//    {
//      const ddr_study_record& study = patient.studies_[j];
//      for(k=0; k < study.series_.size();k++)
//      {
//        const ddr_series_record& series = study.series_[k];
//        if(modlist.find(series.modality_) != string::npos) {
//          isPatientValid = true;
//          break;
//        }
//      }
//    }
//    if(!isPatientValid)
//      continue;
//
//    // Добавляем запись о пациенте
//    toCatalog.push_back(ddr_patient_record());
//    ddr_patient_record& patientNew = toCatalog.back();
//    patientNew.id_ = patient.id_;
//    patientNew.name_ = patient.name_;
//    patientNew.birthDate_ = patient.birthDate_;
//    patientNew.sex_ = patient.sex_;
//
//    // Проверка наличия нужных серий в исследовании
//    bool isStudyValid = false;
//    for(j=0; j < patient.studies_.size(); j++)
//    {
//      const ddr_study_record& study = patient.studies_[j];
//      for(k=0; k < study.series_.size();k++)
//      {
//        const ddr_series_record& series = study.series_[k];
//        if(modlist.find(series.modality_) != string::npos) {
//          isStudyValid = true;
//          break;
//        }
//      }
//      if(!isStudyValid)
//        continue;
//
//      // Добавляем запись об исследовании
//      patientNew.studies_.push_back(ddr_study_record());
//      ddr_study_record& studyNew = patientNew.studies_.back();
//      studyNew.uid_ = study.uid_;
//      studyNew.id_ = study.id_;
//      studyNew.institution_ = study.institution_;
//      studyNew.datetime_ = study.datetime_;
//      studyNew.number_ = study.number_;
//      studyNew.descr_ = study.descr_;
//      studyNew.patient_ = &patientNew;
//
//      // Добавляем серии
//      for(k=0; k < study.series_.size();k++)
//      {
//        const ddr_series_record& series = study.series_[k];
//        if(modlist.find(series.modality_) != string::npos)
//          studyNew.series_.push_back(series);
//      }
//    }
//  }
//}
