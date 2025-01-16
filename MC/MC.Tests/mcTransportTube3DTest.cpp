// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcVRMLDumper.h"
#include "../MC/mcTransportTube3D.h"
#include <memory>
#include <fstream>

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
			p.u = geomVector3D(0, 0.9, -1.45) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
			// 1.3
			p.u = geomVector3D(5, 0, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.4
			p.u = geomVector3D(0, 0, -4) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());

			// 2.1
			p.p = geomVector3D(0.5, 1.5, -1);
			p.u = geomVector3D(0, 0.85, -1.75) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
			// 2.2
			p.u = geomVector3D(0, 0.75, -0.5) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
			// 2.3
			p.u = geomVector3D(0, 1, 1.25) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());

			// 3.1
			p.p = geomVector3D(-0.2, 2, 3);
			p.u = geomVector3D(0, 1.45, 2.25) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
			// 3.2
			p.u = geomVector3D(0, 1.45, 2.5) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
			// 3.3
			p.u = geomVector3D(0, 0, 2.5) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());
		}

		TEST_METHOD(getDistanceInside)
		{
			auto errmsg = L"getDistanceInside failed";
			double d, dexp = 0;
			auto transport = createTestTransport();
			mcParticle p;

			// 1.1 (входной торец под углом)
			p.p = geomVector3D(-0.3, 0.5, -1);
			p.u = geomVector3D(0.3, -0.5, -2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.2 (боковая стенка)
			p.u = geomVector3D(0, 0.75, -0.75) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			//Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::AreEqual(dexp, d, 0.1, errmsg, LINE_INFO());
			// 1.3 (боковой сегмент с переходом в сегмент вперед)
			p.u = geomVector3D(0, 0, 1.6757) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, 0.001, errmsg, LINE_INFO());

			// 2.1
			p.p = geomVector3D(0.2, 1, 1.5);
			p.u = geomVector3D(0.0, 1, 2.45) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
			// 2.2
			p.u = geomVector3D(0, -0.3, 1.5) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());

			// 3.1
			p.p = geomVector3D(-0.2, 1.25, 2);
			p.u = geomVector3D(0.0, 1.45, 2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
			// 3.2
			p.u = geomVector3D(0, 1.25, 1.5) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInside(p);
			Assert::AreEqual(dexp, d, 0.2, errmsg, LINE_INFO());
		}

	private:
		static std::unique_ptr<mcTransportTube3D> createTestTransport()
		{
			std::vector<geomVector3D> pts;
			std::vector<double> rs;

			pts.push_back(geomVector3D(0, 1.5, 2.5));
			rs.push_back(0.0);

			pts.push_back(geomVector3D(0, 1, 2));
			rs.push_back(0.5);

			pts.push_back(geomVector3D(0, 0, 1));
			rs.push_back(0.5);

			pts.push_back(geomVector3D(0, 0, -2));
			rs.push_back(1.0);

			auto t = std::make_unique<mcTransportTube3D>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), pts, rs);
			t->setColor(0.8, 0.8, 0.8, 0.0);

			std::ofstream os("c:/tmp/Tube3D.wrl");
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
