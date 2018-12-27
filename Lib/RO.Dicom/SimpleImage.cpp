#include "pch.h"
#include "SimpleImage.h"
#include "DicomObject.h"
#include "DataSet.h"
#include "DicomDictionary.h"
#include "ElementFactory.h"
//#include "Element.h"
#include "FileMetaInformation.h"
#include "Tag.h"
//#include "TransferSyntax.h"
#include "Uid.h"
#include <fstream>
//#include <strstream>

using namespace RO::Dicom;

std::string CreateUid(int n)
{
	const char* uid_root = "1.2.826.0.1.3680043.8.326.1";

	SYSTEMTIME t;
	GetLocalTime(&t);

	unsigned long tstamp =
		((unsigned long(t.wYear - 2000) << 26) & 0xFC000000) +
		((unsigned long(t.wMonth) << 22) & 0x03C00000) +
		((unsigned long(t.wDay) << 17) & 0x003E0000) +
		((unsigned long(t.wHour) << 12) & 0x0001F000) +
		((unsigned long(t.wMinute) << 6) & 0x00000FC0) +
		(unsigned long(t.wSecond) & 0x0000003F);

	char s[128];
	sprintf_s(s, "%s.%i.%i.%i.%li.%i", uid_root, 0, 0, 999, tstamp, n);
	return s;
}

SimpleImage::SimpleImage(void) : 
nx_(0), ny_(0), 
psx_(0), psy_(0), 
x0_(0), y0_(0), z0_(0), 
pixels_(nullptr)
{
}

SimpleImage::~SimpleImage(void)
{
	if (pixels_ != nullptr)
		delete [] pixels_;
}

void SimpleImage::SetModality(const char* modality)
{
	modality_ = modality;
}

void SimpleImage::SetPixelData(unsigned short* data, double rs, double ri)
{
	if (pixels_ != nullptr)
		delete [] pixels_;
	pixels_ = new unsigned short[nx_ * ny_];
	if (rs == 1.0 && ri == 0)
		memcpy(pixels_, data, nx_ * ny_ * sizeof(unsigned short));
	else
	{
		unsigned short *pdest = pixels_;
		for (int j = 0; j < ny_; j++)
		{
			short* mm = (short*)pixels_ + nx_ * j;
			for (int i = 0; i < nx_; i++)
			{
				double pixel = double(*mm++)*rs + ri;
				if (pixel<0)
					pixel = 0;
				*pdest++ = (unsigned short) pixel;
			}
		}
	}
}

void FillFromDcm(SimpleImage* obj, DicomObject& dcm)
{
	// Наполнение объекта изображения
	auto dset = dcm.Dset;

	obj->SetModality(dset->GetElementString(Tag::Modality).c_str());
	obj->SetResolution(dset->GetElement(Tag::Columns)->GetInt()[0], dset->GetElement(Tag::Rows)->GetInt()[0]);

	auto el_ps = dset->GetElement(Tag::ImagePlanePixelSpacing);
	if (el_ps.get() != nullptr)
	{
		auto f = el_ps->GetDouble();
		if (f.size() == 2)
			obj->SetPixelSize(f[0], f[1]);
	}

	auto el_pos = dset->GetElement(Tag::RTImagePosition);
	if (el_pos.get() != nullptr)
	{
		auto f = el_pos->GetDouble();
		if (f.size() == 2)
			obj->SetImageCorner(f[0] - obj->Psx() * 0.5, f[1] - obj->Psy() * 0.5);
	}

	//// Положение среза по оси Z
	//void SetPosition(double z) { z0_ = z; }

	std::string modality = dcm.Dset->GetElement(Tag::Modality)->GetString();

	auto ri_elem = dcm.Dset->GetElement(Tag::RescaleIntercept);
	double ri = ri_elem == nullptr ? 0 : ri_elem->GetDouble()[0];

	auto rs_elem = dcm.Dset->GetElement(Tag::RescaleSlope);
	double rs = rs_elem == nullptr ? 1.0 : rs_elem->GetDouble()[0];

	if (modality == "CT")
		ri += 1024;

	void* data = dcm.PixelData();
	if (data != nullptr)
		obj->SetPixelData((unsigned short*) data, rs, ri);
}

void SimpleImage::ReadFromStream(std::istream& is)
{
	DicomObject dcm;
	dcm.ReadFromStream(is);
	FillFromDcm(this, dcm);
}

void SimpleImage::ReadFromFile(const char* fname)
{
	DicomObject dcm;
	dcm.ReadFromFile(fname);
	FillFromDcm(this, dcm);
}

void SimpleImage::WtiteToStream(std::ostream& os) const
{
	// Наполнение объекта DICOM из объекта изображения.
	DicomObject dcm;
	auto dset = std::make_shared<DataSet>();

	// (M) Patient
	dset->AddElementString(Tag::PatientName, "RR Portal Image");
	dset->AddElementString(Tag::PatientID, "1111");
	//	obj->addElement(&dcmElement(DCM_PATBIRTHDATE, rtuW32DateConverter(patient_->birthday_).getDateDicomString().c_str()));
	//	obj->addElement(&dcmElement(DCM_PATSEX, patient_->sex_.c_str()));
	//	obj->addElement(&dcmElement(DCM_PATOTHERNAMES, patient_->givenname_.c_str()));
	//	obj->addElement(&dcmElement(DCM_PATCOMMENTS, patient_->comment_.c_str()));

	// (M) General Study
	SYSTEMTIME t;
	GetLocalTime(&t);
	dset->AddElementDate(Tag::StudyDate, t.wYear, t.wMonth, t.wDay);
	dset->AddElementTime(Tag::StudyTime, t.wHour, t.wMinute, t.wSecond, t.wMilliseconds);

	// (M) RT Series
	dset->AddElementString(Tag::Modality, modality_.c_str());
	dset->AddElementString(Tag::StudyInstanceUID, CreateUid(10).c_str());
	dset->AddElementString(Tag::SeriesInstanceUID, CreateUid(20).c_str());
	dset->AddElementString(Tag::SeriesNumber, "1");
	dset->AddElementString(Tag::SeriesDescription, "Simple image");
	dset->AddElementString(Tag::PatientPosition, "HFS");

	// (M) General Equipment
	dset->AddElementString(Tag::Manufacturer, "ROISS");
	dset->AddElementString(Tag::ManufacturerModelName, "RO.Dicom");
	dset->AddElementString(Tag::SoftwareVersions, "V.1.0");

	// (M) General Image
	dset->AddElementString(Tag::InstanceNumber, "1");
	dset->AddElementString(Tag::PatientOrientation, "");
	dset->AddElementString(Tag::ImageType, "DERIVED\\SECONDARY\\DRR");
	dset->AddElementDate(Tag::AcquisitionDate, t.wYear, t.wMonth, t.wDay);
	dset->AddElementTime(Tag::AcquisitionTime, t.wHour, t.wMinute, t.wSecond, t.wMilliseconds);

	// (M) Image Pixel
	dset->AddElementInt(Tag::SamplesPerPixel, 1);
	dset->AddElementString(Tag::PhotometricInterpretation, "MONOCHROME2");
	dset->AddElementInt(Tag::Rows, ny_);
	dset->AddElementInt(Tag::Columns, nx_);
	dset->AddElementInt(Tag::BitsAllocated, 16);
	dset->AddElementInt(Tag::BitsStored, 16);
	dset->AddElementInt(Tag::HighBit, 15);
	dset->AddElementInt(Tag::PixelRepresentation, 0);
	dset->AddElementString(Tag::WindowCenter, "16000");
	dset->AddElementString(Tag::WindowWidth, "32000");

	auto ditem = DicomDictionary::Item(Tag::PixelData);
	auto v = std::make_shared<ElementOtherString>(Tag::PixelData, ditem->Vr);
	v->SetData((const char*) pixels_, nx_ * ny_ * sizeof(unsigned short));
	dset->AddElement(v);

	// (M) RT Image
	std::vector<double> vd2(2);
	vd2[0] = psx_; vd2[1] = psy_;
	dset->AddElementDouble(Tag::ImagePlanePixelSpacing, vd2);

	vd2[0] = x0_ + psx_ * 0.5; vd2[1] = y0_ + psy_ * 0.5;
	dset->AddElementDouble(Tag::RTImagePosition, vd2);

	// (M) SOP Common
	dset->AddElementString(Tag::SOPClassUID, RTImageStorage);
	dset->AddElementString(Tag::SOPInstanceUID, CreateUid(30).c_str());
	dset->AddElementDate(Tag::InstanceCreationDate, t.wYear, t.wMonth, t.wDay);
	dset->AddElementTime(Tag::InstanceCreationTime, t.wHour, t.wMinute, t.wSecond, t.wMilliseconds);
	dset->AddElementString(Tag::SpecificCharacterSet, "ISO_IR 100");

	// Запись объекта DICOM
	dcm.MetaInfo = FileMetaInformation::CreateDefault();
	dcm.Dset = dset;
	dcm.WtiteToStream(os);
}

void SimpleImage::WtiteToFile(const char* fname) const
{
	std::ofstream os(fname, std::ios::binary);
	if (os.fail())
		throw std::exception((std::string("Can not create file: \"") + fname + "\"").c_str());
	this->WtiteToStream(os);
	os.close();
}
