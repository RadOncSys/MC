#include "pch.h"
#include "ElementFactory.h"

using namespace RO::Dicom;
using namespace RO::Dicom::Vr;

std::shared_ptr<Element> ElementFactory::MakeElement(UINT32 tag, VR vr)
{
	std::shared_ptr<Element> v;

	if(vr == VR::AE)
		v = std::make_shared<ElementString>(tag, vr);
	else if(vr == VR::AS)
		v = std::make_shared<ElementAS>(tag, vr);
	else if(vr == VR::AT)
		v = std::make_shared<ElementAT>(tag, vr);
	else if(vr == VR::CS)
		v = std::make_shared<ElementCS>(tag, vr);
	else if(vr == VR::DA)
		v = std::make_shared<ElementDA>(tag, vr);
	else if(vr == VR::TM)
		v = std::make_shared<ElementTM>(tag, vr);
	else if(vr == VR::DT)
		v = std::make_shared<ElementDT>(tag, vr);
		
	// user will parse text for numerical value himself
	else if(vr == VR::DS || vr == VR::IS || vr == VR::LO || vr == VR::LT
		 || vr == VR::PN || vr == VR::SH || vr == VR::ST )	
		v = std::make_shared<ElementString>(tag, vr);

	// user will parse text for numerical value himself
	else if(vr == VR::OB || vr == VR::OF || vr == VR::OW || vr == VR::UN || vr == VR::UT)	
		v = std::make_shared<ElementOtherString>(tag, vr);

	else if(vr == VR::FL)
		v = std::make_shared<ElementFL>(tag, vr);
	else if(vr == VR::FD)
		v = std::make_shared<ElementFD>(tag, vr);

	else if(vr == VR::SL)
		v = std::make_shared<ElementSL>(tag, vr);
	else if(vr == VR::SS)
		v = std::make_shared<ElementSS>(tag, vr);
	else if(vr == VR::UI)
		v = std::make_shared<ElementUI>(tag, vr);
	else if(vr == VR::UL)
		v = std::make_shared<ElementUL>(tag, vr);
	else if(vr == VR::US)
		v = std::make_shared<ElementUS>(tag, vr);
	else if(vr == VR::SQ)
		v = std::make_shared<ElementSQ>(tag, vr);
	else
		throw std::exception("unknown representatrion");

	return v;
}
