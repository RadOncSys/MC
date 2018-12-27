#include "pch.h"
#include "DicomObject.h"
#include "DataSet.h"
#include "Element.h"
#include "FileMetaInformation.h"
#include "Tag.h"
#include <fstream>
#include <strstream>

using namespace RO::Dicom;

DicomObject::DicomObject(void)
{
}

DicomObject::~DicomObject(void)
{
}

void DicomObject::ReadFromStream(std::istream& is, bool skipPixels)
{
	//MetaInfo = std::make_shared<FileMetaInformation>(*dset.get());
	MetaInfo = std::make_shared<FileMetaInformation>();
	MetaInfo->ReadFromStream(is);
	auto ts = MetaInfo->Ts();

	auto dset = std::make_shared<DataSet>(ts);
	dset->Read(is);

	//MetaInfo->RemoveFromDataSet(*dset.get());
	Dset = dset;
}

void DicomObject::ReadFromFile(std::string fname, bool skipPixels)
{
	std::ifstream is(fname.c_str(), std::ios::binary);
	if(is.fail()) 
		throw std::exception(("File:\"" + fname + "\" not found").c_str());
	
	// Skip preamble
	is.seekg(128);

	char prefix[4];
	is.read(prefix, 4);
	if(std::string(prefix, 4) != "DICM")
		throw std::exception(("File:\"" + fname + "\" is not DICOM file").c_str());

	ReadFromStream(is, skipPixels);
	is.close();
}

void DicomObject::ReadFromBuffer(const char* data, unsigned int len)
{
	std::istrstream is(data, len);
	if(is.fail()) 
		throw std::exception("Can not create stream from memory buffer");
	
	// Skip preamble
	is.seekg(128);

	char prefix[4];
	is.read(prefix, 4);
	if(std::string(prefix, 4) != "DICM")
		throw std::exception("Memory buffer do not contain DICOM data");

	ReadFromStream(is);
}

void DicomObject::WtiteToStream(std::ostream& os)
{
	if (!MetaInfo)
		throw std::exception("Dicom object do not have meta information to write to file");
	std::string s0(128, (char)0);
	os.write(s0.c_str(), s0.length());
	os.write("DICM", 4);

	MetaInfo->WriteToStream(os);
	Dset->Write(os, MetaInfo->Ts());
}

void DicomObject::WtiteToFile(std::string fname)
{
	std::ofstream os(fname.c_str(), std::ios::binary);
	if(os.fail()) 
		throw std::exception(("Can not create file: \"" + fname + "\"").c_str());
	WtiteToStream(os);
	os.close();
}

void* DicomObject::PixelData()
{
	if(Dset.get() == nullptr)
		return nullptr;
	auto pe = std::dynamic_pointer_cast<ElementOtherString>(Dset->GetElement(Tag::PixelData));
	if(pe == nullptr)
		return nullptr;
	const char* data = pe->GetData();
	//return static_cast<void*>(data);
	return (void*)(data);
}

long DicomObject::PixelDataLength()
{
	if(Dset.get() == nullptr)
		return 0;
	auto pe = std::dynamic_pointer_cast<ElementOtherString>(Dset->GetElement(Tag::PixelData));
	if(pe == nullptr)
		return 0;
	return pe->GetLength();
}
