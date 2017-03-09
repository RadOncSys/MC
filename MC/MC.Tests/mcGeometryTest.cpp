#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcGeometry.h"
#include "../geometry/vec3d.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcGeometryTest)
	{
	public:

		TEST_METHOD(getDistanceToInfiniteCylinderInside)
		{
			double d;
			d = mcGeometry::getDistanceToInfiniteCylinderInside(geomVector3D(0, 0, 10), geomVector3D(1, 0, 0), 10);
			Assert::AreEqual(10, d, DBL_EPSILON, L"getDistanceToInfiniteCylinderInside failed", LINE_INFO());

			d = mcGeometry::getDistanceToInfiniteCylinderInside(geomVector3D(5, 0, -100), geomVector3D(-1, 0, 0), 10);
			Assert::AreEqual(15, d, DBL_EPSILON, L"getDistanceToInfiniteCylinderInside failed", LINE_INFO());

			d = mcGeometry::getDistanceToInfiniteCylinderInside(geomVector3D(5, 0, -100), geomVector3D(0, -1, 0), 10);
			Assert::AreEqual(sqrt(100.0 - 25.0), d, DBL_EPSILON, L"getDistanceToInfiniteCylinderInside failed", LINE_INFO());

			d = mcGeometry::getDistanceToInfiniteCylinderInside(geomVector3D(5, 5, -100), geomVector3D(0, 0, -1), 10);
			Assert::AreEqual(DBL_MAX, d, DBL_EPSILON, L"getDistanceToInfiniteCylinderInside failed", LINE_INFO());

			geomVector3D v(0, -1, -1);
			v.normalize();
			d = mcGeometry::getDistanceToInfiniteCylinderInside(geomVector3D(5, 0, -100), v, 10);
			Assert::AreEqual(sqrt((100.0 - 25.0) * 2), d, DBL_EPSILON * 10, L"getDistanceToInfiniteCylinderInside failed", LINE_INFO());
		}

		TEST_METHOD(getDistanceToSphereInside)
		{
			double r = 10.0, d;
			geomVector3D p(0, 0, 0);
			geomVector3D v(1, 0, 0);

			d = mcGeometry::getDistanceToSphereInside(p, v, r);
			Assert::AreEqual(r, d, DBL_EPSILON, L"getDistanceToSphereInside failed", LINE_INFO());

			p(1) = -r * 0.5;
			d = mcGeometry::getDistanceToSphereInside(p, v, r);
			Assert::AreEqual(r * sqrt(0.75), d, 10 * DBL_EPSILON, L"getDistanceToSphereInside failed", LINE_INFO());

			p(1) = -r * 0.5;
			p(2) = -r * 0.5;
			d = mcGeometry::getDistanceToSphereInside(p, v, r);
			Assert::AreEqual(r * sqrt(0.5), d, 10 * DBL_EPSILON, L"getDistanceToSphereInside failed", LINE_INFO());
		}

		TEST_METHOD(getDistanceToSphereOutside)
		{
			double r = 10.0, d;
			geomVector3D p(-2*r, 0, 0);
			geomVector3D v(1, 0, 0);

			p(1) = -r * 0.5;
			d = mcGeometry::getDistanceToSphereOutside(p, v, r);
			Assert::AreEqual(r * (2.0 - sqrt(0.75)), d, 10 * DBL_EPSILON, L"getDistanceToSphereOutside failed", LINE_INFO());

			p(1) = -r * 0.5;
			p(2) = -r * 0.5;
			d = mcGeometry::getDistanceToSphereOutside(p, v, r);
			Assert::AreEqual(r * (2.0 - sqrt(0.5)), d, 10 * DBL_EPSILON, L"getDistanceToSphereOutside failed", LINE_INFO());
		}
	};
}
