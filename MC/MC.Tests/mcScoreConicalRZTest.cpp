#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcScoreConicalRZ.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcScoreConicalRZTest)
	{
	public:

		TEST_METHOD(ScorePoint)
		{
			auto score = createTestScore();
			// 1
			score->ScorePoint(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(0, 0, 10.1));
			Assert::AreEqual(1.0, score->Dose(0, 50), DBL_EPSILON * 10, L"ScorePoint failed", LINE_INFO());
			// 2
			score->ScorePoint(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(0, 2.0, 5.1));
			Assert::AreEqual(1.0, score->Dose(21, 25), DBL_EPSILON * 10, L"ScorePoint failed", LINE_INFO());
			// 3
			score->ScorePoint(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(5.0, 0, 20.1));
			Assert::AreEqual(1.0, score->Dose(44, 100), DBL_EPSILON * 10, L"ScorePoint failed", LINE_INFO());
		}

		TEST_METHOD(ScoreLine)
		{
			int i, j;
			auto score = createTestScore();
			// 1 - вдоль центральной оси на длину скоринга
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(0.05, 0, 0), geomVector3D(0.05, 0, 35.0));
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(0.05, 0, 35.0), geomVector3D(0.05, 0, 0));
			for (i = 0; i < 175; i++)
				Assert::AreEqual(2.0 / 175, score->Dose(0, i), DBL_EPSILON * 10, L"ScoreLine failed", LINE_INFO());
			
			// 2 - вдоль радиуса между границей и центром в двух напралениях
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(0, 0, 20.1), geomVector3D(0, 5.0*90.1 / 80.0, 20.1));
			score->ScoreLine(1.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(0, 5.0*90.1 / 80.0, 20.1), geomVector3D(0, 0, 20.1));
			Assert::AreEqual(2.0 / 50 + 2.0 / 175, score->Dose(0, 100), DBL_EPSILON * 100, L"ScoreLine failed", LINE_INFO());
			for (i = 1; i < 50; i++)
				Assert::AreEqual(2.0 / 50, score->Dose(i, 100), DBL_EPSILON * 100, L"ScoreLine failed", LINE_INFO());
			
			// 3 - по ассиметричной диагонали в двух направлениях
			score->ScoreLine(2.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(5.0 * 70.0 / 80.0, 0, 0), geomVector3D(0, 5.0 * 105.0 / 80.0, 35.0));
			score->ScoreLine(2.0, 0, *(mcRegionReference*)nullptr, mc_particle_t::MCP_PHOTON, geomVector3D(5.0 * 105.0 / 80.0, 0, 35.0), geomVector3D(0, 5.0 * 70.0 / 80.0, 0));

			double sum = 0;
			for (i = 0; i < 50; i++)
			{
				for (j = 0; j < 175; j++)
					sum += score->Dose(i, j);
			}

			Assert::AreEqual(4.0 / 175, score->Dose(49, 174), DBL_EPSILON * 1000, L"ScoreLine failed", LINE_INFO());
			Assert::AreEqual(8.0, sum, DBL_EPSILON * 1000, L"ScoreLine failed", LINE_INFO());
		}

	private:
		static std::shared_ptr<mcScoreConicalRZ> createTestScore()
		{
			return std::shared_ptr<mcScoreConicalRZ>(
				new mcScoreConicalRZ("mcScoreConicalRZ test", 1, 50, 175, 5.0, 0.0, 35.0, 10.0, 80.0));
		}
	};
}
