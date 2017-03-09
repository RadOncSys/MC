#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcSourceLEBA.h"
#include "../MC/mcThread.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcSourceLEBATest)
	{
	private:
		int nThreads_ = 1;
		double emean_ = 6.0;
		double ewidth_ = 2.0;
		double z_ = 0.0;
		double rsigma_ = 0.1;

	public:
		TEST_METHOD(sample)
		{
			double se = ewidth_ / sqrt(12.0);
			EBATypeTest(L"PRISM distributed electron beam", spectrum_distr_t::SPECTRUM_PRISM, se);
		
			//se = ewidth_ * sqrt(5.0 / 72.0);
			se = ewidth_ / sqrt(24.0);
			EBATypeTest(L"TRIANGLE distributed electron beam", spectrum_distr_t::SPECTRUM_TRIANGLE, se);

			se = ewidth_ / 2.0;
			EBATypeTest(L"GAUSS distributed electron beam", spectrum_distr_t::SPECTRUM_GAUSS, se);
		}

	private:
		void EBATypeTest(const wchar_t* name, spectrum_distr_t sptype, double se_expected)
		{
			// Expected sigmar
			double sigmar_expected = rsigma_ * sqrt(6.0);

			mcSourceLEBA source("Test LEBA source", nThreads_, sptype, profile_distr_t::PROFILE_EXPONENT, emean_, ewidth_, z_, rsigma_);

			mcThread thread;
			thread.setId(0);
			mcParticle p;

			double ecount = 0;
			double e2count = 0;
			double r2count = 0;

			int i, n = 100000;
			for (i = 0; i < n; i++)
			{
				source.sample(p, &thread);
				ecount += p.ke;
				e2count += (p.ke - emean_) * (p.ke - emean_);
				r2count += p.p.sqLengthXY();
			}

			// Средняя энергия при всех распределениях должна быть равна средней заданной в описании спектра.
			Assert::AreEqual(ecount / n, emean_, 5.0 / sqrt(double(n)), name, LINE_INFO());
			Assert::AreEqual(sqrt(e2count / n), se_expected, se_expected * 5.0 / sqrt(double(n)), name, LINE_INFO());

			// Проверка самплинга положения
			Assert::AreEqual(sqrt(r2count / n), sigmar_expected, sigmar_expected * 5.0 / sqrt(double(n)), name, LINE_INFO());
		}
	};
}
