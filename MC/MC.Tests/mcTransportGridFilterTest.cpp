// Radiation Oncology Monte Carlo open source project
//
// Author: [2024] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcVRMLDumper.h"
#include "../MC/mcTransportGridFilter.h"
#include <memory>
#include <fstream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportGridFilterTest)
	{
	public:

		TEST_METHOD(getDistanceOutside)
		{
			auto errmsg = L"getDistanceOutside failed";
			double d = 0, dexp = 0;
			auto transport = createTestTransport();
			mcParticle p;

			double h = transport->az();
			double a = transport->ax() / 2;

			//
			// Тестирование расстояния снаружи - это расстояние до параллелпипеда
			// 

			// 1.1
			p.p = geomVector3D(-a / 3, -a / 3, 2 * h);
			p.u = geomVector3D(-a / 2, -a / 2, h) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 1.2
			p.u = geomVector3D(a * 2, -a / 2, h) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, 0.2, errmsg, LINE_INFO());
			// 1.3
			p.u = geomVector3D(-a / 2, -a / 2, 2 * h) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());

			// 2.1
			p.p = geomVector3D(-1.5 * a, -a / 3, 0.75 * h);
			p.u = geomVector3D(-a, -a / 2, 0.25 * h) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 2.2
			p.u = geomVector3D(-1.1 * a, -a / 2, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, 0.2, errmsg, LINE_INFO());
			// 2.3
			p.u = geomVector3D(-2 * a, -a / 2, 0.25 * h) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());

			// 3.1
			p.p = geomVector3D(-a / 2, -a / 3, -h);
			p.u = geomVector3D(-a / 3, -a / 2, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			// 3.2
			p.u = geomVector3D(a * 3, -a / 2, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, 0.2, errmsg, LINE_INFO());
			// 3.3
			p.u = geomVector3D(-a / 3, -a / 2, -2 * h) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceOutside(p);
			Assert::AreEqual(DBL_MAX, d, TEST_EPSILON, errmsg, LINE_INFO());

		}

		TEST_METHOD(getIdxAtPoint)
		{
			auto errmsg = L"getIdxAtPoint failed";
			double d, dexp = 0;
			auto transport = createTestTransport();
			mcParticle p;
			bool isInBrick;
			int idx;

			double x0 = transport->psx() / 2;
			double y0 = transport->psy() / 2;
			double psz = transport->psz();

			int nx = transport->nx();
			int ny = transport->ny();
			int nz = transport->nz();

			int i = nx / 2;
			int j = ny / 2;
			int k = 0;

			// Тестируем 3 брикета (два крайних и по середине
			double bs1 = transport->getBrickSize(0);
			double bs2 = transport->getBrickSize(nz / 3);
			double bs3 = transport->getBrickSize(nz - 1);

			double z1 = psz / 3;
			double z2 = psz * (nz / 3 + 0.5);
			double z3 = psz * (nz - 0.3);

			// 1
			k = 0;
			p.p = geomVector3D(x0 - bs1 * 1.0001 / 2, y0 / 2, z1);
			idx = transport->getIdxAtPoint(p.p, p.region.gidx_, isInBrick);
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			Assert::IsFalse(isInBrick, errmsg, LINE_INFO());
			// 2
			p.p = geomVector3D(x0 - bs1 * 0.99 / 2, y0 / 2, z1);
			idx = transport->getIdxAtPoint(p.p, p.region.gidx_, isInBrick);
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			Assert::IsTrue(isInBrick, errmsg, LINE_INFO());

			// 3
			k = nz / 3;
			p.p = geomVector3D(x0 - bs2 * 1.1 / 2, y0, z2);
			idx = transport->getIdxAtPoint(p.p, p.region.gidx_, isInBrick);
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			Assert::IsFalse(isInBrick, errmsg, LINE_INFO());
			// 4
			p.p = geomVector3D(x0, y0 + bs2 * 0.5 / 2, z2);
			idx = transport->getIdxAtPoint(p.p, p.region.gidx_, isInBrick);
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			Assert::IsTrue(isInBrick, errmsg, LINE_INFO());

			// 5
			k = nz - 1;
			p.p = geomVector3D(x0, y0 - bs3 * 1.1 / 2, z3);
			idx = transport->getIdxAtPoint(p.p, p.region.gidx_, isInBrick);
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			Assert::IsFalse(isInBrick, errmsg, LINE_INFO());
			// 6
			p.p = geomVector3D(x0 - bs3 * 0.7 / 2, y0, z3);
			idx = transport->getIdxAtPoint(p.p, p.region.gidx_, isInBrick);
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			Assert::IsTrue(isInBrick, errmsg, LINE_INFO());
		}

		TEST_METHOD(getDistanceInsideVoxel)
		{
			auto errmsg = L"getDistanceInsideVoxel failed";
			double d, dexp = 0;
			auto transport = createTestTransport();
			mcParticle p;
			bool isHitCell;
			int idx;
			short gidx[3];

			double x0 = transport->psx() / 2;
			double y0 = transport->psy() / 2;
			double psz = transport->psz();

			int nx = transport->nx();
			int ny = transport->ny();
			int nz = transport->nz();

			int i = nx / 2;
			int j = ny / 2;
			int k = 0;

			// Тестируем 3 брикета (два крайних и по середине
			double bs1 = transport->getBrickSize(0);
			double bs2 = transport->getBrickSize(nz / 3);
			double bs3 = transport->getBrickSize(nz - 1);

			double z1 = psz / 3;
			double z2 = psz * (nz / 3 + 0.5);
			double z3 = psz * (nz - 0.3);

			// 1.1 - строго вниз
			k = 0;
			p.region.gidx_[0] = i; p.region.gidx_[1] = j; p.region.gidx_[2] = k;
			p.region.idx_ = (k * ny + j) * nx + i;
			p.region.subidx_ = 0;

			p.p = geomVector3D(x0 - bs1 * 1.0001 / 2, y0 / 2, z1);
			p.u = geomVector3D(x0 - bs1 * 1.0001 / 2, y0 / 2, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, -1, errmsg, LINE_INFO());
			// 1.2 - строго вверх
			p.u = geomVector3D(x0 - bs1 * 1.0001 / 2, y0, psz) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k + 1) * ny + j) * nx + i, errmsg, LINE_INFO());
			// 1.3 - в сторону брикета
			p.u = geomVector3D(x0 - bs1 / 2, y0, z1 / 2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsFalse(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			// 1.4 - от брикета
			p.u = geomVector3D(0, y0, psz * 0.9) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i - 1, errmsg, LINE_INFO());

			// 2.1 - вниз
			p.region.gidx_[0] = i; p.region.gidx_[1] = j; p.region.gidx_[2] = k;
			p.region.idx_ = (k * ny + j) * nx + i;
			p.region.subidx_ = 1;

			p.p = geomVector3D(x0 - bs1 * 0.99 / 2, y0 / 2, z1);
			p.u = geomVector3D(x0 - bs1 / 2, y0, 0) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, -1, errmsg, LINE_INFO());
			// 2.2 - вверх
			p.u = geomVector3D(x0, y0 - bs1 / 2, psz) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k + 1) * ny + j) * nx + i, errmsg, LINE_INFO());
			// 2.3 - в сторону полости
			p.u = geomVector3D(x0 + bs1 / 2, y0 / 2, z1 / 2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsFalse(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i + 1, errmsg, LINE_INFO());

			// 3.1 - строго вниз
			k = nz / 3;
			p.region.gidx_[0] = i; p.region.gidx_[1] = j; p.region.gidx_[2] = k;
			p.region.idx_ = (k * ny + j) * nx + i;
			p.region.subidx_ = 0;

			p.p = geomVector3D(x0 - bs2 * 1.1 / 2, y0 / 2, z2);
			p.u = geomVector3D(x0 - bs2 * 1.1 / 2, y0 / 2, psz * k) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k - 1) * ny + j) * nx + i, errmsg, LINE_INFO());
			// 3.2 - строго вверх
			p.u = geomVector3D(x0 - bs2 * 1.1 / 2, y0, psz * (k + 1)) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k + 1) * ny + j) * nx + i, errmsg, LINE_INFO());
			// 3.3 - в сторону брикета
			p.u = geomVector3D(x0 - bs2 / 2, y0, z2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsFalse(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			// 3.4 - от брикета
			p.u = geomVector3D(0, y0 / 2, psz * (k + 0.9)) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i - 1, errmsg, LINE_INFO());

			// 4.1 - вниз
			p.region.gidx_[0] = i; p.region.gidx_[1] = j; p.region.gidx_[2] = k;
			p.region.idx_ = (k * ny + j) * nx + i;
			p.region.subidx_ = 1;

			p.p = geomVector3D(x0 - bs2 / 4, y0 / 2, z2);
			p.u = geomVector3D(x0 - bs2 / 2, y0, psz * k) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k - 1) * ny + j) * nx + i, errmsg, LINE_INFO());
			// 4.2 - вверх
			p.u = geomVector3D(x0, y0 - bs2 / 4, psz * (k + 1)) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k + 1) * ny + j) * nx + i, errmsg, LINE_INFO());
			// 4.3 - в сторону полости
			p.u = geomVector3D(x0 + bs2 / 2, y0 / 2, z2) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsFalse(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i + 1, errmsg, LINE_INFO());

			// 5.1 - строго вниз
			k = nz - 1;
			p.region.gidx_[0] = i; p.region.gidx_[1] = j; p.region.gidx_[2] = k;
			p.region.idx_ = (k * ny + j) * nx + i;
			p.region.subidx_ = 0;

			p.p = geomVector3D(x0 - bs3 * 1.1 / 2, y0 / 2, z3);
			p.u = geomVector3D(x0 - bs3 * 1.1 / 2, y0 / 2, psz * k) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k - 1)* ny + j)* nx + i, errmsg, LINE_INFO());
			// 5.2 - строго вверх
			p.u = geomVector3D(x0 - bs3 * 1.1 / 2, y0, psz * (k + 1)) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, -1, errmsg, LINE_INFO());
			// 5.3 - в сторону брикета
			p.u = geomVector3D(x0 - bs3 / 2, y0, z3) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsFalse(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i, errmsg, LINE_INFO());
			// 5.4 - от брикета
			p.u = geomVector3D(0, y0, z3) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i - 1, errmsg, LINE_INFO());

			// 6.1 - вниз
			p.region.gidx_[0] = i; p.region.gidx_[1] = j; p.region.gidx_[2] = k;
			p.region.idx_ = (k * ny + j) * nx + i;
			p.region.subidx_ = 1;

			p.p = geomVector3D(x0 - bs3 * 0.99 / 2, y0 / 2, z3);
			p.u = geomVector3D(x0 - bs3 / 2, y0, psz * k) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, ((k - 1)* ny + j)* nx + i, errmsg, LINE_INFO());
			// 6.2 - вверх
			p.u = geomVector3D(x0, y0 - bs3 / 2, psz * (k + 1)) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsTrue(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, -1, errmsg, LINE_INFO());
			// 6.3 - в сторону полости
			p.u = geomVector3D(x0 + bs3 / 2, y0 / 2, z3) - p.p;
			dexp = p.u.length();
			p.u.normalize();
			d = transport->getDistanceInsideVoxel(p, gidx, idx, isHitCell);
			Assert::AreEqual(dexp, d, TEST_EPSILON, errmsg, LINE_INFO());
			Assert::IsFalse(isHitCell, errmsg, LINE_INFO());
			Assert::AreEqual<int>(idx, (k * ny + j) * nx + i + 1, errmsg, LINE_INFO());
		}

	private:
		static std::shared_ptr<mcTransportGridFilter> createTestTransport()
		{
			int nx = 14;
			int ny = 14;
			int nz = 13;
			double psx = 0.5;
			double psy = 0.5;
			double psz = 0.096;

			auto t = std::make_shared<mcTransportGridFilter>(
				geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), nx, ny, nz, psx, psy, psz);

			t->setBrickSize(0, 0.499, 0.499);
			t->setBrickSize(1, 0.405, 0.405);
			t->setBrickSize(2, 0.375, 0.375);
			t->setBrickSize(3, 0.345, 0.345);
			t->setBrickSize(4, 0.315, 0.315);
			t->setBrickSize(5, 0.290, 0.290);
			t->setBrickSize(6, 0.265, 0.265);
			t->setBrickSize(7, 0.240, 0.240);
			t->setBrickSize(8, 0.220, 0.220);
			t->setBrickSize(9, 0.200, 0.200);
			t->setBrickSize(10, 0.175, 0.175);
			t->setBrickSize(11, 0.140, 0.140);
			t->setBrickSize(12, 0.105, 0.105);

			t->setColor(0.8, 0.8, 0.8, 0.0);

			std::ofstream os("c:/tmp/GridFilter.wrl");
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
