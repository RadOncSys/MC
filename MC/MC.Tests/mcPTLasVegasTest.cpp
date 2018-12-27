// Radiation Oncology Monte Carlo open source project
//
// Author: [2015-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcPTLasVegas.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcPTLasVegasTest)
	{
	public:

		TEST_METHOD(getDistanceOutside)
		{
			auto errmsg = L"getDistanceOutside failed";
			double d;
			mcPTLasVegas transport(geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0));

			// Координаты первого и вторго узлов сетки и параметры первого цилиндра.
			double x0 = transport.Xs()[0], x1 = transport.Xs()[1], y0 = transport.Ys()[0], y1 = transport.Ys()[1];
			double d0 = transport.Ds()[0], h0 = transport.Hs()[0];
			mcParticle p;

			// 1
			p.p = geomVector3D(-2, 5, 10 + transport.Z());
			p.u = geomVector3D(0, 0, -1);
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(10, d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 2
			p.p = geomVector3D(transport.A()/2 + 1, 0, 0);
			p.u = geomVector3D(-sqrt(0.5), 0, sqrt(0.5));
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(sqrt(2.0), d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 3
			p.p = geomVector3D(-transport.A()/2 - 1, 0, 0);
			p.u = geomVector3D(sqrt(0.5), 0, sqrt(0.5));
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(sqrt(2.0), d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 4
			p.p = geomVector3D((x0 + x1)*0.5, (y0 + y1)*0.5, -1);
			p.u = geomVector3D(0, 0, 1);
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(1, d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 5
			p.p = geomVector3D(x0 - d0 - h0, y0, -d0/2);
			p.u = geomVector3D(sqrt(0.5), 0, sqrt(0.5));
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(d0*sqrt(0.5), d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 6
			p.p = geomVector3D(x0, y0, 0);
			p.u = geomVector3D(0, 0, 1);
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(h0, d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 7
			p.p = geomVector3D(x0 + d0 / 2 - 1.5*h0, y0, -h0);
			p.u = geomVector3D(sqrt(0.5), 0, sqrt(0.5));
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(1.5*h0*sqrt(2.0), d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 8
			p.p = geomVector3D(x0, y0, h0);
			p.u = geomVector3D(0, 0, 1);
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(0, d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 9
			p.p = geomVector3D(x0, y0, h0);
			p.u = geomVector3D(0, 0, -1);
			d = transport.getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, DBL_EPSILON * 10, errmsg, LINE_INFO());
		}

		TEST_METHOD(getDistanceInside)
		{
			auto errmsg = L"getDistanceInside failed";
			double d;
			mcPTLasVegas transport(geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0));

			// Координаты первого и вторго узлов сетки и параметры первого цилиндра.
			double x0 = transport.Xs()[0], x1 = transport.Xs()[1], y0 = transport.Ys()[0], y1 = transport.Ys()[1];
			double d0 = transport.Ds()[0], h0 = transport.Hs()[0];
			mcParticle p;

			// 1
			p.p = geomVector3D(0, 0, 1);
			p.u = geomVector3D(0, sqrt(0.5), sqrt(0.5));
			d = transport.getDistanceInside(p);
			Assert::AreEqual((transport.Z()-1)*sqrt(2.0), d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 2
			p.p = geomVector3D(-transport.A()/2-1, 0, 1);
			p.u = geomVector3D(-sqrt(0.5), 0, -sqrt(0.5));
			d = transport.getDistanceInside(p);
			Assert::AreEqual(sqrt(2.0), d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 3
			p.p = geomVector3D(-transport.A() / 2 - 1, 0, 2);
			p.u = geomVector3D(-sqrt(0.5), 0, -sqrt(0.5));
			d = transport.getDistanceInside(p);
			Assert::AreEqual(sqrt(2.0), d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 4
			p.p = geomVector3D(x0 + d0/4 + h0, y0, 2*h0);
			p.u = geomVector3D(-sqrt(0.5), 0, -sqrt(0.5));
			d = transport.getDistanceInside(p);
			Assert::AreEqual(sqrt(2.0)*h0, d, DBL_EPSILON * 10, errmsg, LINE_INFO());
			// 5
			p.p = geomVector3D(x0 + 0.5*d0 + 0.5*h0, y0, h0);
			p.u = geomVector3D(-sqrt(0.5), 0, -sqrt(0.5));
			d = transport.getDistanceInside(p);
			Assert::AreEqual(sqrt(0.5)*h0, d, DBL_EPSILON * 10, errmsg, LINE_INFO());
		}
	};
}
