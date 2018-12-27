// Radiation Oncology Monte Carlo open source project
//
// Author: [2015-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcScoreSphereMatrix.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcScoreSphereMatrixTest)
	{
	public:

		TEST_METHOD(ScorePoint)
		{
			auto score = createTestScore();
			// 1
			score->ScorePoint(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, geomVector3D(0.0001, 1, 100));
			Assert::AreEqual(1.0, score->Dose(0, 8), DBL_EPSILON * 10, L"ScorePoint failed", LINE_INFO());
			// 2
			score->ScorePoint(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, geomVector3D(-0.0001, -1, -100));
			Assert::AreEqual(1.0, score->Dose(17, 26), DBL_EPSILON * 10, L"ScorePoint failed", LINE_INFO());
			// 3
			score->ScorePoint(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, geomVector3D(100.0, 1, 1));
			Assert::AreEqual(1.0, score->Dose(8, 0), DBL_EPSILON * 10, L"ScorePoint failed", LINE_INFO());
		}

		TEST_METHOD(ScoreLine)
		{
			int i, j;
			auto score = createTestScore();
			// 1 - строго по радиусу внутри детектора
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, geomVector3D(0.0001, 1, 99.5), geomVector3D(0.0001, 1, 100.5));
			Assert::AreEqual(1.0, score->Dose(0, 8), DBL_EPSILON * 10, L"ScoreLine failed", LINE_INFO());
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, geomVector3D(0.0001, 1, 100.5), geomVector3D(0.0001, 1, 99.5));
			Assert::AreEqual(2.0, score->Dose(0, 8), DBL_EPSILON * 10, L"ScoreLine failed", LINE_INFO());
			
			// 2 - аналогичный трек строго по радиусу но длиннее с высживанием только половины
			geomVector3D v(0.0001, 1, 99.5);
			v.normalize();
			score->ScoreLine(2.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v * 98.0, v * 102.0);
			Assert::AreEqual(3.0, score->Dose(0, 8), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());
			score->ScoreLine(2.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v * 102.0, v * 98.0);
			Assert::AreEqual(4.0, score->Dose(0, 8), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());

			// 3 - длинный трек со множеством пересечений из пересечений осей X и Z с внешней сферой
			geomVector3D v1(1, 0.0001, 100);
			geomVector3D v2(100, 0.0001, 1);
			v1.normalize();
			v2.normalize();
			score->ScoreLine(100.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v1 * 102, v2 * 102);
			double d00 = 2.0444829574519723;
			Assert::AreEqual(d00, score->Dose(0, 0), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(d00, score->Dose(8, 0), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());
			
			score->ScoreLine(100.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v2 * 102, v1 * 102);
			Assert::AreEqual(2*d00, score->Dose(0, 0), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(2*d00, score->Dose(8, 0), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());

			// 4 - длинный трек со множеством пересечений из пересечений осей X и Y с внешней сферой
			v1 = geomVector3D(1, 100, -0.0001);
			v2 = geomVector3D(100, 1, -0.0001);
			v1.normalize();
			v2.normalize();
			score->ScoreLine(100.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v1 * 102, v2 * 102);
			d00 = 2.0444829574519723;
			d00 = score->Dose(9, 0);
			Assert::AreEqual(d00, score->Dose(9, 0), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(d00, score->Dose(9, 8), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());
			
			score->ScoreLine(100.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v2 * 102, v1 * 102);
			Assert::AreEqual(2*d00, score->Dose(9, 0), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(2*d00, score->Dose(9, 8), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());

			// 5 - Пересечение вокселей паралелей 
			v1 = geomVector3D(100, -0.0001, -5);
			v2 = geomVector3D(100, -0.0001, 5);
			v1.normalize();
			v2.normalize();
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v1 * 101, v2 * 101);
			Assert::AreEqual(0.5, score->Dose(8, 35), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(0.5, score->Dose(9, 35), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());

			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v2 * 101, v1 * 101);
			Assert::AreEqual(1.0, score->Dose(8, 35), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(1.0, score->Dose(9, 35), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());

			// 6 - Пересечение вокселей меридианов 
			v1 = geomVector3D(-100, -5, -0.0001);
			v2 = geomVector3D(-100, 5, -0.0001);
			v1.normalize();
			v2.normalize();
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v1 * 101, v2 * 101);
			Assert::AreEqual(0.5, score->Dose(9, 17), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(0.5, score->Dose(9, 18), DBL_EPSILON * 20, L"ScoreLine failed", LINE_INFO());

			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_NEGATRON, v2 * 101, v1 * 101);
			Assert::AreEqual(1.0, score->Dose(9, 17), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(1.0, score->Dose(9, 18), DBL_EPSILON * 40, L"ScoreLine failed", LINE_INFO());
		}

	private:
		static std::shared_ptr<mcScoreSphereMatrix> createTestScore()
		{
			return std::shared_ptr<mcScoreSphereMatrix>(
				new mcScoreSphereMatrix("mcScoreSphereMatrix test", 1, 18, 36, 99, 101));
		}
	};
}
