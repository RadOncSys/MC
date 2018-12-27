// Radiation Oncology Monte Carlo open source project
//
// Author: [2015-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcTransportCone.h"
#include "../geometry/vec3d.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportConeTest)
	{
	public:

		TEST_METHOD(getDistanceInside)
		{
			// GG 20171126
			// <size unit = "cm" radius = "3.4" height = "1.0" / >
			// Wrong particle transport inside object : Flattenig Filter
			//	Position : 2.10727 - 0.651952       0.351231        1
			//	Direction : 0.748649 - 0.662605 - 0.0219135      1
			//	Dnear : 0
			//	Distance = 1.79769e+308

			mcTransportCone transport(geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), 3.4, 1.0);

			mcParticle p;
			p.p.set(2.10727, -0.651952, 0.351231);
			p.u.set(0.748649, -0.662605, -0.0219135);

			double dist = transport.getDistanceInside(p);

			//Assert::IsTrue(dist > 0, L"getDistanceInside failed", LINE_INFO());
			Assert::AreEqual(0.0, dist, 0.0001, L"getDistanceInside failed", LINE_INFO());
		}
	};
}
