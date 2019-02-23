// Radiation Oncology Monte Carlo open source project
//
// Author: [2015-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcTransportRectanglePolygonSideHole.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportRectanglePolygonSideHoleTest)
	{
	public:
		
		TEST_METHOD(getDistanceOutside)
		{
			double d, dexpected;
			auto transport = createTestTransport();
			mcParticle p;

			// 1.1
			p.p = geomVector3D(-1, 4, -10);
			p.u = geomVector3D(0, 0, 1);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 1.2
			p.u = geomVector3D(-1, 1, 1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(10 * sqrt(3.0), d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 1.3
			p.u = geomVector3D(-2, 0, 1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());

			// 2.1
			p.p = geomVector3D(2, 1, 15);
			p.u = geomVector3D(-4, 0, -11);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 2.1.1
			p.u = geomVector3D(1, 0, -11);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 2.2
			p.u = geomVector3D(1, 0, -1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual((5*sqrt(2.0) + sqrt(0.5)), d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 2.3
			p.u = geomVector3D(10, 0, -5);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 2.4
			p.u = geomVector3D(20, 0, -5);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 2.5
			p.u = geomVector3D(0, 4, -11);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 2.6
			p.u = geomVector3D(0, -5, -13);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());

			// 3.1
			p.p = geomVector3D(1, -1, 6);
			p.u = geomVector3D(-1, 0, 0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(4.0, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 3.2
			p.u = geomVector3D(1, 0, -1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(2 * sqrt(2.0), d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 3.3
			p.u = geomVector3D(0, 0, -1);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 3.4
			p.u = geomVector3D(0, 0, 1);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());

			// 4.1
			p.p = geomVector3D(15, 4, -5);
			p.u = geomVector3D(-2, 0, 7);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
			// 4.2
			p.u = geomVector3D(-17, 0, 6);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());

			// 5.1
			p.p = geomVector3D(-15, -13, 15);
			p.u = geomVector3D(3, 0, -7);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceOutside failed", LINE_INFO());
		}
		
		TEST_METHOD(getDistanceInside)
		{
			double d, dexpected;
			auto transport = createTestTransport();
			mcParticle p;

			// 1.1
			p.p = geomVector3D(-8, -2, 6);
			p.u = geomVector3D(1, 0, 1);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0 * sqrt(0.5), d, TEST_EPSILON, L"getDistanceInside failed", LINE_INFO());
			// 1.2
			p.u = geomVector3D(-3, 0, 4);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceInside failed", LINE_INFO());
			// 1.3
			p.u = geomVector3D(-4, 0, -2);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceInside failed", LINE_INFO());
			// 1.4
			p.u = geomVector3D(4, 0, -6);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceInside failed", LINE_INFO());
			// 1.5
			p.u = geomVector3D(6, 0, -2);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceInside failed", LINE_INFO());

			// 2.1
			p.p = geomVector3D(8, 7, 4);
			p.u = geomVector3D(-5, -2, 0);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceInside failed", LINE_INFO());
			// 2.2
			p.u = geomVector3D(0, 2, 5);
			dexpected = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexpected, d, TEST_EPSILON, L"getDistanceInside failed", LINE_INFO());
		}

	private:
        static std::unique_ptr<mcTransportRectanglePolygonSideHole> createTestTransport()
		{
			double d = 10, x1 = -2, x2 = 3, y1 = -4, y2 = 5;

			std::vector<double> poly_z;
			std::vector<double> poly_x;
			std::vector<double> poly_y;

			poly_z.push_back(0);
			poly_x.push_back(0);
			poly_y.push_back(0);

			poly_z.push_back(5);
			poly_x.push_back(0);
			poly_y.push_back(0);

			poly_z.push_back(10);
			poly_x.push_back(5);
			poly_y.push_back(5);

			auto t = std::make_unique<mcTransportRectanglePolygonSideHole>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0),
				d, d, poly_z, poly_x, poly_y);
			t->SetFieldSize(x1, x2, y1, y2);

			return std::move(t);
		}
	};
}