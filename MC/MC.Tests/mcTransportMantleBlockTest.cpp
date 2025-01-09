// Radiation Oncology Monte Carlo open source project
//
// Author: [2025] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcVRMLDumper.h"
#include "../MC/mcTransportMantleBlock.h"
#include <memory>
#include <fstream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportMantleBlockTest)
	{
	public:

		TEST_METHOD(getDistanceOutside)
		{
			auto errmsg = L"mcTransportTube3DTest::getDistanceOutside failed";
			double d = 0, dexp = 0;
			auto transport = createTestTransport();
			mcParticle p;

			// 1.1 (в обратную сторону)
			p.p = geomVector3D(-0.5, -0.5, 1.75);
			p.u = geomVector3D(0, 1.5, 2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.2 (в отверстие и внутреннюю стенку)
			p.u = geomVector3D(0.5, -0.5, 0.2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.3 (Сквозь отверстие)
			p.u = geomVector3D(0.25, +0.1, 0.0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.4 (в тело)
			p.u = geomVector3D(1, 0.5, 1.0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.5 (мимо объекта)
			p.u = geomVector3D(1.5, 0.5, 1.0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());

			// 2.1 (в обратную сторону)
			p.p = geomVector3D(1, 0.5, -0.5);
			p.u = geomVector3D(0, 1.5, -1) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 2.2 (мимо объекта)
			p.u = geomVector3D(1.5, 0.5, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 2.3 (в тело)
			p.u = geomVector3D(1.25, -0.5, 0.0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());

			// 3.1 (вылет из отверстия вперед)
			p.p = geomVector3D(-0.5, -0.25, 0.375);
			p.u = geomVector3D(-0.75, 0, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 3.2 (вылет из отверстия назад)
			p.u = geomVector3D(0.25, 0.25, 1.0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 3.3 (в боковую стенку)
			p.u = geomVector3D(-1, -0.75, 0.1) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());

			// 4.1 (мимо объекта)
			p.p = geomVector3D(-2.75, -0.5, 0.375);
			p.u = geomVector3D(-1.5, 0.5, 1) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 4.2 (мимо объекта)
			p.u = geomVector3D(-1.5, -1, 0.0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 4.3 (в боковую стенку)
			p.u = geomVector3D(-1.5, 0, 0.9) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
		}

		TEST_METHOD(getDistanceInside)
		{
			auto errmsg = L"mcTransportTube3DTest::getDistanceInside failed";
			double d, dexp = 0;
			auto transport = createTestTransport();
			mcParticle p;

			// 1.1 (в боковую стенку цилиндра)
			p.p = geomVector3D(0.75, 0.75, 0.625);
			p.u = geomVector3D(0, 1.5, 0.375) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.2 (в нижний торец)
			p.u = geomVector3D(-0.5, 0.5, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.3 (во внутреннее отверстие)
			p.u = geomVector3D(-0.5, 0.3, 0.1) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.4 (во внутреннее отверстие)
			p.u = geomVector3D(0.5, -0.3, 0.2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
		}

	private:
		static std::unique_ptr<mcTransportMantleBlock> createTestTransport()
		{
			std::vector<double> px;
			std::vector<double> py;

			px.push_back(-1);
			py.push_back(-0.8);

			px.push_back(-1);
			py.push_back(0.3);

			px.push_back(0.5);
			py.push_back(0.3);

			px.push_back(0.5);
			py.push_back(-0.8);

			auto t = std::make_unique<mcTransportMantleBlock>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), 1.5, 1.0, px, py);
			t->setColor(1.0, 0.5, 0, 0.1);

			std::ofstream os("c:/tmp/MantleBlock.wrl");
			if (!os.fail())
			{
				mcVRMLDumper::dumpHead(os);
				mcVRMLDumper::dumpWorldAxis(os);
				t->dumpVRML(os);
			}
			return std::move(t);
		}
	};
}
