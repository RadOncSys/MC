// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcTransportJawPairRounded.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportJawPairRoundedTest)
	{
	public:

		TEST_METHOD(getDistanceOutside)
		{
			double d;
			auto transport = createTestTransport();
			mcParticle p;

			// 1
			p.p = geomVector3D(-30, 0, 5);
			p.u = geomVector3D(1, 0, 0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(2.0, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 2
			p.p = geomVector3D(20, 0, 5);
			p.u = geomVector3D(-1, 0, 0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(3.0, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());

			// 3
			p.p = geomVector3D(-14, 20, 5);
			p.u = geomVector3D(0, -1, 0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(5.1, d, 0.2, L"getDistanceOutside failed", LINE_INFO());
			// 4
			p.p = geomVector3D(3, -20, 5);
			p.u = geomVector3D(0, 1, 0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(5.1, d, 0.2, L"getDistanceOutside failed", LINE_INFO());

			// 5
			p.p = geomVector3D(-3, 0, 5);
			p.u = geomVector3D(-1, -1, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(5 * sqrt(2.0), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 6
			p.p = geomVector3D(-8, 0, 5);
			p.u = geomVector3D(1, 1, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(5 * sqrt(2.0), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());

			// 6
			p.p = geomVector3D(2, 0, 11);
			p.u = geomVector3D(-1, 0, -1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(sqrt(2.0), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 7
			p.p = geomVector3D(17, 0, -1);
			p.u = geomVector3D(-1, 0, 1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(sqrt(2.0), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
		}

		TEST_METHOD(getDistanceInside)
		{
			double d;
			auto transport = createTestTransport();
			mcParticle p;

			// 1 (1)
			p.p = geomVector3D(-13, 0, 5);
			p.u = geomVector3D(-1, 1, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(15.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 2 (1)
			p.u = geomVector3D(1, 1, -1);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5 * sqrt(3.0), d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 1 (2)
			p.p = geomVector3D(2, 0, 5);
			p.u = geomVector3D(-1, 1, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5 * sqrt(2.0), d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 2 (2)
			p.u = geomVector3D(1, -1, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(15.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 2 (3)
			p.u = geomVector3D(0, -1, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(15.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 2 (4)
			p.u = geomVector3D(1, 1, 1);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5 * sqrt(3.0), d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());
		}

	private:
		static std::shared_ptr<mcTransportJawPairRounded> createTestTransport()
		{
			double radius = 15.0;
			double height = 10.0;
			double width = 5.0;
			double fsx1 = -8.0;
			double fsx2 = -3.0;

			auto t = std::make_shared<mcTransportJawPairRounded>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), radius, height, width);
			t->setFS(fsx1, fsx2);

			return t;
		}
	};
}