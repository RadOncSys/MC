// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcETransportConvexPolygonCircle.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcETransportConvexPolygonCircleTest)
	{
	public:

		TEST_METHOD(getDistanceOutside)
		{
			double d;
			auto transport = createTestTransport();
			mcParticle p;
			
			// 5 (3)
			p.p = geomVector3D(-2, 0, -4);
			p.u = geomVector3D(2, 0, 4);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(1.0/p.u.z(), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 6 (3)
			p.p = geomVector3D(0, -2, -4);
			p.u = geomVector3D(0, -2, 1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			
			// 7 (4)
			p.p = geomVector3D(-4, 0, 7);
			p.u = geomVector3D(4, 0, -3);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(-2.0 / p.u.z(), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			// 8 (4)
			p.u = geomVector3D(0.5, 0, -3);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(2.8, d, 0.2, L"getDistanceOutside failed", LINE_INFO());
			// 9 (4)
			p.u = geomVector3D(-1, 0, -3);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
			
			// 10 (5)
			p.p = geomVector3D(0, -7, 4.5);
			p.u = geomVector3D(0, 7, -0.5);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(3.5, d, 0.2, L"getDistanceOutside failed", LINE_INFO());
			// 11 (5)
			p.u = geomVector3D(0, 7, -4.5);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(3.0 / p.u.y(), d, DBL_EPSILON * 10, L"getDistanceOutside failed", LINE_INFO());
		}

		TEST_METHOD(getDistanceInside)
		{
			double d;
			auto transport = createTestTransport();
			mcParticle p;

			// 1 (1)
			p.p = geomVector3D(0, 0, 1);
			p.u = geomVector3D(0, 0, -1);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(4.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 2 (1)
			p.u = geomVector3D(0, 2.000000001, 4);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(sqrt(20), d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 3 (2)
			p.p = geomVector3D(2, 0, 1);
			p.u = geomVector3D(1, 0, 0);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(2.0, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());

			// 4 (2)
			p.u = geomVector3D(2, 0, 4);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(3.5777087639996634, d, DBL_EPSILON * 10, L"getDistanceInside failed", LINE_INFO());


			// Ðואכüםûי מבתוךע ס נואכüםûלט ןנמבכולאלט
			auto t2 = createTestTransport2();
			p.p = geomVector3D(-6.3251, 6.01896, -5.84373);
			p.u = geomVector3D(-0.801706, 0.536857, -0.262776);
			d = t2->getDistanceInside(p);
			Assert::AreEqual(0.0, d, 0.1, L"getDistanceInside failed", LINE_INFO());
		}

	private:
		static std::shared_ptr<mcETransportConvexPolygonCircle> createTestTransport()
		{
			std::vector<double> poly_z;
			std::vector<double> poly_r;

			poly_z.push_back(-3); poly_r.push_back(3);
			poly_z.push_back(-1); poly_r.push_back(4);
			poly_z.push_back(4); poly_r.push_back(4);
			poly_z.push_back(5); poly_r.push_back(2);

			return std::make_shared<mcETransportConvexPolygonCircle>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0),
				poly_z, poly_r);
		}

		static std::shared_ptr<mcETransportConvexPolygonCircle> createTestTransport2()
		{
			std::vector<double> poly_z;
			std::vector<double> poly_r;

			poly_z.push_back(-22); poly_r.push_back(5.5);
			poly_z.push_back(-2); poly_r.push_back(9.5);
			poly_z.push_back(5.3); poly_r.push_back(9.5);
			poly_z.push_back(9.8); poly_r.push_back(5);

			return std::make_shared<mcETransportConvexPolygonCircle>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0),
				poly_z, poly_r);
		}
	};
}