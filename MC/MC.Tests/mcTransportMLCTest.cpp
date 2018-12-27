// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcTransportMLC.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportMLCTest)
	{
	public:

		TEST_METHOD(getDistanceOutside)
		{
			double d;
			auto transport = createTestTransport();
			mcParticle p;
			std::wstring err(L"mcTransportMLCTest::getDistanceOutside failed");

			// 1 (0)
			p.p = geomVector3D(0, 0, 5);
			p.u = geomVector3D(1, 0, 0);
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(5.0, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 2 (0)
			p.u = geomVector3D(1, -2, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(2.5*sqrt(5.0)*0.9, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			p.u = geomVector3D(2, 1, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(2.5 * sqrt(5.0), d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 3 (0)
			p.p = geomVector3D(0, 0, 10.0);
			p.u = geomVector3D(10.0*(1-sqrt(1-0.0625)), 7.5, -2.5);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(7.9, d, 0.1, err.c_str(), LINE_INFO());

			// 4 (0)
			p.p = geomVector3D(0, 0, 5);
			p.u = geomVector3D(1, 0, 0.5);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(6.0, d, 0.1, err.c_str(), LINE_INFO());

			p.u = geomVector3D(-1, 0, -0.5);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(6.0, d, 0.1, err.c_str(), LINE_INFO());

			// 5 (0)
			p.u = geomVector3D(0, 0, 1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, 0, err.c_str(), LINE_INFO());
			
			p.u = geomVector3D(0, 0, -1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, 0, err.c_str(), LINE_INFO());

			// 6 (1)
			p.p = geomVector3D(36, 11, 5);
			p.u = geomVector3D(-1, -11, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(11.05, d, 0.05, err.c_str(), LINE_INFO());

			// 7 (1)
			p.u = geomVector3D(-3, -6, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(3.25*sqrt(5), d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 8 (1)
			p.u = geomVector3D(-6, -3, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(3.0*sqrt(5), d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 9 (1)
			p.u = geomVector3D(-20, -2, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(sqrt(402.25), d, 0.05, err.c_str(), LINE_INFO());
			
			p.u = geomVector3D(-40, -2, 0);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(sqrt(1602.25), d, 0.05, err.c_str(), LINE_INFO());

			// 10 (2)
			p.p = geomVector3D(36, 0, 11);
			p.u = geomVector3D(-3, 0, -1);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(sqrt(10.0), d, 0.05, err.c_str(), LINE_INFO());

			// 11 (2)
			p.u = geomVector3D(-1, 0, -6);
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(sqrt(37.0), d, 0.05, err.c_str(), LINE_INFO());

			//
			// Версия с ассиметричным полем
			//
			auto transportAssym = createTestTransportAssym();

			// 12 (3)
			p.p = geomVector3D(-40, -11, 5);
			p.u = geomVector3D(10, 3, 0);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(sqrt(109), d, 0.05, err.c_str(), LINE_INFO());

			// 13 (3)
			p.u = geomVector3D(15, 10, 0);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(sqrt(325), d, 0.05, err.c_str(), LINE_INFO());

			// 14 (3)
			p.u = geomVector3D(10, 15, 0);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(15.5/p.u.y(), d, 0.05, err.c_str(), LINE_INFO());

			p.u = geomVector3D(10, 17, 0);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(sqrt(100+17*17), d, 0.05, err.c_str(), LINE_INFO());

			// 15 (4)
			p.p = geomVector3D(10, 0, 8);
			p.u = geomVector3D(0, 1, 0);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(5*0.84, d, 0.05, err.c_str(), LINE_INFO());

			p.u = geomVector3D(0, -1, 0);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(5 * 0.84, d, 0.05, err.c_str(), LINE_INFO());

			p.u = geomVector3D(0, 0, 1);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, 0, err.c_str(), LINE_INFO());

			p.u = geomVector3D(0, 0, -1);
			p.u.normalize();
			d = transportAssym->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, 0, err.c_str(), LINE_INFO());
		}

		TEST_METHOD(getDistanceInside)
		{
			double d;
			auto transport = createTestTransport();
			mcParticle p;
			std::wstring err(L"mcTransportMLCTest::getDistanceInside failed");

			// 1 (0)
			p.p = geomVector3D(-30, 0, 5);
			p.u = geomVector3D(-1, 0, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 2 (0)
			p.u = geomVector3D(1, 0, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(25.0, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 3 (0)
			p.u = geomVector3D(-1, -1, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0 * 0.9 * sqrt(2), d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			p.u = geomVector3D(-1, 1, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0 * 0.9 * sqrt(2), d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 4 (0)
			p.u = geomVector3D(1, 10, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(sqrt(101) * 0.9, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());
			
			p.u = geomVector3D(-1, 10, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(sqrt(101) * 0.45, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 5 (1)
			p.p = geomVector3D(10, 0, 5);
			p.u = geomVector3D(-1, 0, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 6 (1)
			p.u = geomVector3D(1, 0, 0);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(25.0, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			// 7 (1)
			p.u = geomVector3D(-1, 0, 1);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.8, d, .1, err.c_str(), LINE_INFO());
			
			p.u = geomVector3D(0, 0, -1);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(5.0, d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());

			p.u = geomVector3D(1, 1, 1);
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(sqrt(75), d, DBL_EPSILON * 10, err.c_str(), LINE_INFO());
		}

	private:
		static std::shared_ptr<mcTransportMLC> createTestTransport()
		{
			double focus = 50.0;
			double radius = 10.0;
			double height = 10.0;
			double width = 20.0;
			double length = 30.0;
			double x1 = -5.0;
			double x2 = 5.0;
			double y1 = -5.0;
			double y2 = 5.0;

			auto t = std::make_shared<mcTransportMLC>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), focus, radius, height, width, length, x1, x2, y1, y2);

			return t;
		}

		static std::shared_ptr<mcTransportMLC> createTestTransportAssym()
		{
			double focus = 50.0;
			double radius = 10.0;
			double height = 10.0;
			double width = 20.0;
			double length = 30.0;
			double x1 = 5.0;
			double x2 = 15.0;
			double y1 = -5.0;
			double y2 = 5.0;

			auto t = std::make_shared<mcTransportMLC>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), focus, radius, height, width, length, x1, x2, y1, y2);

			return t;
		}
	};
}