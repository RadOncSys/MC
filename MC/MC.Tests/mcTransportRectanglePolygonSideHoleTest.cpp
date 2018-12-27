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
			double d;
			auto transport = createTestTransport();
			mcParticle p;
			// 1
			p.p = geomVector3D(-2,7,-10);
			p.u = geomVector3D(0,0,1);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, DBL_EPSILON*10, L"getDistanceOutside failed", LINE_INFO());
			// 2
			p.p = geomVector3D(0,12,-10);
			p.u = geomVector3D(0,0,1);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(10.0, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 3
			p.p = geomVector3D(-7,0,15);
			p.u = geomVector3D(0,0,-1);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(8.0, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 4
			p.p = geomVector3D(2,5,5);
			p.u = geomVector3D(0,0,1);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 5
			p.p = geomVector3D(-20,0,3);
			p.u = geomVector3D(1,0,0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(5, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 6
			p.p = geomVector3D(2,0,2);
			p.u = geomVector3D(-1,0,0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(7.0, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 7
			p.p = geomVector3D(3,0,7);
			p.u = geomVector3D(1,0,0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(4.0, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 8
			p.p = geomVector3D(7,0,-10);
			p.u.set(-1, 0, 1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(12.0*sqrt(2.0), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 9
			p.p = geomVector3D(-17,28,1);
			p.u.set(1.0, -0.2, 0.3);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
		}
		
		TEST_METHOD(getDistanceInside)
		{
			double d;
			auto transport = createTestTransport();
			mcParticle p;

			p.p = geomVector3D(7,7,0);
			p.u = geomVector3D(0,0,1);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(7.0, d, DBL_EPSILON*10, L"getDistanceInside failed", LINE_INFO());

			p.p = geomVector3D(7,7,5);
			p.u = geomVector3D(0,0,-1);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			p.p = geomVector3D(7,7,6);
			p.u = geomVector3D(-1,0,0);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(1.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			p.p = geomVector3D(7,7,3);
			p.u = geomVector3D(1,0,0);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(8.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			p.p = geomVector3D(7,17,7);
			p.u = geomVector3D(0,-1,0);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			p.p = geomVector3D(12,0,2);
			p.u.set(-1, 0, 1);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0*sqrt(2.0), d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			p.p = geomVector3D(-10,0,5);
			p.u.set(20, 0, -5);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(sqrt(25.0 + 1.25*1.25), d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());
		}

	private:
        static std::shared_ptr<mcTransportRectanglePolygonSideHole> createTestTransport()
		{
			double d = 10, cx = 5, cy = 10;
			std::vector<double> poly_z;
			std::vector<double> poly_x;
			std::vector<double> poly_y;

			poly_z.push_back(0);
			poly_x.push_back(0+cx);
			poly_y.push_back(0+cy);

			poly_z.push_back(5);
			poly_x.push_back(0+cx);
			poly_y.push_back(0+cy);

			poly_z.push_back(10);
			poly_x.push_back(5+cx);
			poly_y.push_back(5+cy);

			return std::make_shared<mcTransportRectanglePolygonSideHole>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0),
				d + cx, d + cy, poly_z, poly_x, poly_y);
		}
	};
}