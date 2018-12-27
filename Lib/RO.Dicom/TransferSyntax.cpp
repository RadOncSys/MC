#include "pch.h"
#include "TransferSyntax.h"
#include "Uid.h"
#include <map>

using namespace RO::Dicom;
using namespace std;

static std::map<std::string, std::shared_ptr<TransferSyntaxItem>> tsyntaxes_;

shared_ptr<TransferSyntaxItem> RO::Dicom::TransferSyntax::Item(const char* uid)
{
	Populate();

	// Function map::at() throws exception if not found key, wich is not what we expect
	auto it = tsyntaxes_.find(uid);
	if(it == tsyntaxes_.end())
		throw exception("unknown transfer syntax");
	else
		return it->second;
}

void RO::Dicom::TransferSyntax::Populate()
{
	if(tsyntaxes_.empty())
	{
		tsyntaxes_.insert(make_pair(ImplicitVRLittleEndian, std::make_shared<TransferSyntaxItem>("Implicit VR Little Endian: Default Transfer Syntax for DICOM", ImplicitVRLittleEndian, true, false, false, false, false, false)));

		tsyntaxes_.insert(make_pair(ExplicitVRLittleEndian, std::make_shared<TransferSyntaxItem>("Explicit VR Little Endian", ExplicitVRLittleEndian, true, false, true, false, false, false)));

		tsyntaxes_.insert(make_pair(DeflatedExplicitVRLittleEndian, std::make_shared<TransferSyntaxItem>("Deflated Explicit VR Little Endian", DeflatedExplicitVRLittleEndian, true, false, true, true, false, true)));

		tsyntaxes_.insert(make_pair(ExplicitVRBigEndian, std::make_shared<TransferSyntaxItem>("Explicit VR Big Endian", ExplicitVRBigEndian, false, false, true, false, false, false)));

		tsyntaxes_.insert(make_pair(JPEGBaseline1, std::make_shared<TransferSyntaxItem>("JPEG Baseline (Process 1): Default Transfer Syntax for Lossy JPEG 8 Bit Image Compression",
		                                                                                JPEGBaseline1, true, true, true, false, true, false)));

		tsyntaxes_.insert(make_pair(JPEGExtended24, std::make_shared<TransferSyntaxItem>("JPEG Extended (Process 2 &amp; 4): Default Transfer Syntax for Lossy JPEG 12 Bit Image Compression (Process 4 only)",
		                                                                                 JPEGExtended24, true, true, true, false, true, false)));

		tsyntaxes_.insert(make_pair(JPEGLosslessNonHierarchical14, std::make_shared<TransferSyntaxItem>("JPEG Lossless, Non-Hierarchical (Process 14)", JPEGLosslessNonHierarchical14, true, true, true, false, false, true)));

		tsyntaxes_.insert(make_pair(JPEGLossless, std::make_shared<TransferSyntaxItem>("JPEG Lossless, Non-Hierarchical, First-Order Prediction (Process 14 [Selection Value 1]): Default Transfer Syntax for Lossless JPEG Image Compression",
		                                                                               JPEGLossless, true, true, true, false, false, true)));

		tsyntaxes_.insert(make_pair(JPEGLSLossless, std::make_shared<TransferSyntaxItem>("JPEG-LS Lossless Image Compression", JPEGLSLossless, true, true, true, false, false, true)));

		tsyntaxes_.insert(make_pair(JPEGLSLossyNearLossless, std::make_shared<TransferSyntaxItem>("JPEG-LS Lossy (Near-Lossless) Image Compression", JPEGLSLossyNearLossless, true, true, true, false, true, false)));

		tsyntaxes_.insert(make_pair(JPEG2000LosslessOnly, std::make_shared<TransferSyntaxItem>("JPEG 2000 Image Compression (Lossless Only)", JPEG2000LosslessOnly, true, true, true, false, false, true)));

		tsyntaxes_.insert(make_pair(JPEG2000, std::make_shared<TransferSyntaxItem>("JPEG 2000 Image Compression", JPEG2000, true, true, true, false, true, false)));

		tsyntaxes_.insert(make_pair(JPEG2000Part2MultiComponentLosslessOnly, std::make_shared<TransferSyntaxItem>("JPEG 2000 Part 2 Multi-component  Image Compression (Lossless Only)",
		                                                                                                          JPEG2000Part2MultiComponentLosslessOnly, true, true, true, false, false, true)));

		tsyntaxes_.insert(make_pair(JPEG2000Part2MultiComponent, std::make_shared<TransferSyntaxItem>("JPEG 2000 Part 2 Multi-component  Image Compression", JPEG2000Part2MultiComponent, true, true, true, false, true, false)));

		tsyntaxes_.insert(make_pair(JPIPReferenced, std::make_shared<TransferSyntaxItem>("JPIP Referenced", JPIPReferenced, true, false, true, false, false, false)));

		tsyntaxes_.insert(make_pair(JPIPReferencedDeflate, std::make_shared<TransferSyntaxItem>("JPIP Referenced Deflate", JPIPReferencedDeflate, true, false, true, true, false, true)));

		tsyntaxes_.insert(make_pair(MPEG2, std::make_shared<TransferSyntaxItem>("MPEG2 Main Profile @ Main Level", MPEG2, true, true, true, false, true, false)));

		////#define                     "1.2.840.10008.1.2.4.101"
		//		tsyntaxes_.insert(make_pair(MPEG2MainProfileHighLevel, shared_ptr<TransferSyntaxItem>(
		//			new TransferSyntaxItem("", MPEG2MainProfileHighLevel, true, false, true, false, false, false))));

		////#define                "1.2.840.10008.1.2.4.102"
		//		tsyntaxes_.insert(make_pair(MPEG4AVCH264HighProfileLevel41, shared_ptr<TransferSyntaxItem>(
		//			new TransferSyntaxItem("", MPEG4AVCH264HighProfileLevel41, true, false, true, false, false, false))));

		////#define    "1.2.840.10008.1.2.4.103"
		//		tsyntaxes_.insert(make_pair(MPEG4AVCH264BDCompatibleHighProfileLevel41, shared_ptr<TransferSyntaxItem>(
		//			new TransferSyntaxItem("", MPEG4AVCH264BDCompatibleHighProfileLevel41, true, false, true, false, false, false))));

		tsyntaxes_.insert(make_pair(RLELossless, std::make_shared<TransferSyntaxItem>("RLE Lossless", RLELossless, true, true, true, false, false, true)));

		tsyntaxes_.insert(make_pair(RFC2557MIMEEncapsulation, std::make_shared<TransferSyntaxItem>("RFC 2557 MIME encapsulation", RFC2557MIMEEncapsulation, true, false, true, false, false, false)));

		tsyntaxes_.insert(make_pair(XMLEncoding, std::make_shared<TransferSyntaxItem>("XML Encoding", XMLEncoding, true, false, true, false, false, false)));
	}
}

