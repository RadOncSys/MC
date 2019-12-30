#include "pch.h"
#include "FileMetaInformation.h"
#include "Element.h"
#include "ElementFactory.h"
#include "Vr.h"
#include "Tag.h"
#include "TransferSyntax.h"
#include "Uid.h"
#include <istream>
#include <exception>

using namespace RO::Dicom;

FileMetaInformation::FileMetaInformation()
{
}

//FileMetaInformation::FileMetaInformation(DataSet& data)
//{
//	const char* errMsg = "Dataset does not have metafile information";
//
//	data.CopyElementToDataSet(elements_, FileMetaInformationGroupLength, true, errMsg);	// (0002,0000)
//	data.CopyElementToDataSet(elements_, FileMetaInformationVersion, true, errMsg);		// (0002,0001)
//	data.CopyElementToDataSet(elements_, MediaStorageSOPClassUID, true, errMsg);		// (0002,0002)
//	data.CopyElementToDataSet(elements_, MediaStorageSOPInstanceUID, true, errMsg);		// (0002,0003)
//	data.CopyElementToDataSet(elements_, TransferSyntaxUID, true, errMsg);				// (0002,0010)
//	data.CopyElementToDataSet(elements_, ImplementationClassUID, true, errMsg);			// (0002,0012)
//	data.CopyElementToDataSet(elements_, ImplementationVersionName, false, errMsg);		// (0002,0013)
//	data.CopyElementToDataSet(elements_, SourceApplicationEntityTitle, false, errMsg);	// (0002,0016)
//	data.CopyElementToDataSet(elements_, PrivateInformationCreatorUID, false, errMsg);	// (0002,0100)
//	data.CopyElementToDataSet(elements_, PrivateInformation, false, errMsg);			// (0002,0102)
//}

void FileMetaInformation::ReadFromStream(std::istream& is)
{
	const char* errmsg = "error in DataSet parsing from stream";
	auto ts = TransferSyntax::Item(ExplicitVRLittleEndian);

	while(true)
	{
		UINT32 nextValue = Element::PeekUINT32(is);
		UINT16 group = TagGroup(nextValue);
		if(group > 2)
			break;

		UINT32 tag = DataSet::ReadTag(is);
		if(is.fail()) throw std::exception(errmsg);

		Vr::VR vr = DataSet::ReadVr(is);
		if(is.fail()) throw std::exception(errmsg);

		UINT16 len = DataSet::ReadLength(is);
		if(is.fail()) throw std::exception(errmsg);

		//auto pos = is.tellg();

		auto v = ElementFactory::MakeElement(tag, vr);
		v->Read(is, len, ts);

		elements_.AddElement(v);
	}
}

void FileMetaInformation::WriteToStream(std::ostream& os)
{
	auto ts = TransferSyntax::Item(ExplicitVRLittleEndian);

	auto elements = elements_.Elements();
	for each(auto element in *elements)
	{
		DataSet::WriteTag(os, element.first);
		DataSet::WriteVr(os, element.second->vr());
		element.second->Write(os, ts);
	}
}

//void FileMetaInformation::RemoveFromDataSet(DataSet& data)
//{
//	auto elements = elements_.Elements();
//	std::for_each(elements.begin(), elements.end(),
//		[&data](std::pair<UINT32, std::shared_ptr<Element>> x)
//	{
//		data.DeleteElement(x.second->tag());
//	});
//}

std::shared_ptr<TransferSyntaxItem> FileMetaInformation::Ts()
{
	auto element = elements_.GetElement(TransferSyntaxUID);
	if (element)
	{
		std::string uid = element->GetString();
		ts_ = TransferSyntax::Item(uid.c_str());
	}
	else
		ts_ = TransferSyntax::Item(ExplicitVRLittleEndian);
	return ts_;
}

std::shared_ptr<FileMetaInformation> FileMetaInformation::CreateDefault()
{
	auto metaset = std::make_shared<FileMetaInformation>();
	DataSet& metaelements = metaset->Elements();

	//metaelements.AddElementString(MakeTag(0x0002, 0x0000), "");

	//std::vector<int> version(2);
	//version[0] = 0; version[1] = 1;
	//metaelements.AddElementInt(MakeTag(0x0002, 0x0001), version);

	//metaelements.AddElementString(MakeTag(0x0002, 0x0002), "1.2.840.10008.5.1.4.1.1.12.1");
	//metaelements.AddElementString(MakeTag(0x0002, 0x0003), "1.3.12.2.1107.5.4.5.28034.30000012120405135550000000213");

	metaelements.AddElementString(Tag::TransferSyntaxUID, ExplicitVRLittleEndian);
	metaelements.AddElementString(MakeTag(0x0002, 0x0012), "1.3.6.1.4.1.25403.1.1.1");
	metaelements.AddElementString(MakeTag(0x0002, 0x0013), "Dicom 0.1");

	//metaelements.AddElementString(MakeTag(0x0002, 0x0016), "RADARCHIVEDCM");
	//metaelements.AddElementString(MakeTag(0x0002, 0x0100), "");
	//metaelements.AddElementString(MakeTag(0x0002, 0x0102), "");

	return metaset;
}
