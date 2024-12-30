// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcTransportTube3D.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportTube3DTest)
	{
	public:

		TEST_METHOD(getDistanceOutside)
		{
			auto errmsg = L"getDistanceOutside failed";
			double d = 0, dexp = 0;
			auto transport = createTestTransport();
			mcParticle p;

			// 1.1
			p.p = geomVector3D(-0.5, 2, -3);
			p.u = geomVector3D(0.5, 0.5, -2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.2
			p.u = geomVector3D(-0.5, 0.5, -1) - p.p;
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(2.0, d, TEST_EPSILON, errmsg, LINE_INFO());


		}

		TEST_METHOD(getDistanceInside)
		{
			auto errmsg = L"getDistanceInside failed";
			double d;
			auto transport = createTestTransport();
			mcParticle p;

			// 1.1
			p.p = geomVector3D(0, 0, 1);
			p.u = geomVector3D(0, 0, -1) - p.p;
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(3.0, d, TEST_EPSILON, errmsg, LINE_INFO());


		}

	private:
		static std::unique_ptr<mcTransportTube3D> createTestTransport()
		{
			std::vector<geomVector3D> pts;
			std::vector<double> rs;

			pts.push_back(geomVector3D(0, 0, -2));
			rs.push_back(1.0);

			pts.push_back(geomVector3D(0, 0, 1));
			rs.push_back(0.5);

			pts.push_back(geomVector3D(0, 1, 2));
			rs.push_back(0.5);

			pts.push_back(geomVector3D(0, 1.5, 2.5));
			rs.push_back(0);

			auto t = std::make_unique<mcTransportTube3D>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), pts, rs);

			return std::move(t);
		}
	};
}
