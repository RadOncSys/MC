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

			//
			// Дополнительный тест в связи с проблемами при ранспорте в брахитерапевтическом источнике, 
			// воспроизводящий условия проблемного случая
			//
			transport = createTestTransport2();
			p.p = geomVector3D(0.012850686228402975, 0.019238845775936410, 0.17900099635016914);
			p.u = geomVector3D(0.046703267315407700, 0.071450298763663683, 0.99635016918182373);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(0.0384, d, 0.01, errmsg, LINE_INFO());

			p.p = geomVector3D(-0.00900814, 0.000891232, -0.179001);
			p.u = geomVector3D(0.671535, -0.0894564, -0.735553);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(0.0798, d, 0.01, errmsg, LINE_INFO());

			p.p = geomVector3D(-0.015729, -0.000711092, -0.179001);
			p.u = geomVector3D(-0.854681, -0.0385922, -0.517717);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(0.0312, d, 0.01, errmsg, LINE_INFO());

			p.p = geomVector3D(0.00587706, -0.000194463, -0.179001);
			p.u = geomVector3D(-0.326909, 0.0918469, -0.940582);
			d = transport->getDistanceInside(p);
			Assert::AreEqual(0.1493, d, 0.01, errmsg, LINE_INFO());
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

		static std::unique_ptr<mcTransportTube3D> createTestTransport2()
		{
			std::vector<geomVector3D> pts;
			std::vector<double> rs;

			pts.push_back(geomVector3D(0, 0, 0.242));
			rs.push_back(0.0);

			pts.push_back(geomVector3D(0, 0, 0.20));
			rs.push_back(0.045);

			pts.push_back(geomVector3D(0, 0, -2.1));
			rs.push_back(0.045);

			auto t = std::make_unique<mcTransportTube3D>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), pts, rs);
			t->setColor(0.8, 0.8, 0.8, 0.0);

			std::ofstream os("c:/tmp/Tube3D_2.wrl");
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
