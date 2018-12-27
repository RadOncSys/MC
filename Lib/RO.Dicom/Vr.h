// THIS CODE BELONGS TO Radiation Oncology Intellectual Systems and Services LLC
// Copyright (c) 2013, ROISS LLC. All rights reserved
//
// Author: Gennady Gorlachev (ggorlachev@roiss.ru)
//---------------------------------------------------------------------------
#pragma once

namespace RO
{
	namespace Dicom
	{
		namespace Vr
		{
			// PS 3.5-2011 Page 24 
			// Table 6.2-1 DICOM VALUE REPRESENTATIONS

			enum VR {
				AE = 0x4541, //!< Application Entity
				AS = 0x5341, //!< Age String
				AT = 0x5441, //!< Attribute Tag
				CS = 0x5343, //!< Code String
				DA = 0x4144, //!< Date
				DS = 0x5344, //!< Decimal String
				DT = 0x5444, //!< Date Time
				FL = 0x4C46, //!< Floating point single
				FD = 0x4446, //!< Floating point double
				IS = 0x5349, //!< Integer String
				LO = 0x4f4c, //!< Long string
				LT = 0x544c, //!< Long Text
				OB = 0x424f, //!< Other Byte String
				OF = 0x464f, //!< Other Float String
				OW = 0x574f, //!< Other Word String
				PN = 0x4e50, //!< Person Name
				SH = 0x4853, //!< Short String
				SL = 0x4C53, //!< Signed long
				SQ = 0x5153, //!< Sequence
				SS = 0x5353, //!< Signed Short
				ST = 0x5453, //!< Short text
				TM = 0x4d54, //!< Time
				UI = 0x4955, //!< Unique Identifier
				UL = 0x4C55, //!< Unsigned Long
				UN = 0x4e55, //!< Unknown
				US = 0x5355, //!< Unsigned Short
				UT = 0x5455  //!< Unlimited Text
			};
		}
	}
}
