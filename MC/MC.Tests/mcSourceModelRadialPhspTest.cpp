// Radiation Oncology Monte Carlo open source project
//
// Author: [2015-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcSourceModelRadialPhsp.h"
//#include "../MC/mcThread.h"
#include "../MC/mcRng.h"
//#include <memory>
#include <strstream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcSourceModelRadialPhspTest)
	{
	private:
		const int np_ = 20;
		double sample1_ [20] = { 30258, 30581.1, 31213.6, 32132.7, 33378.1, 35072.9, 37129.2, 39590.3, 42499.7, 45890.1, 49554.2, 53541.5, 57637.6, 61798.7, 65372, 68092.5, 69955.9, 71332.8, 72242.8, 72726.2 };
		double sample2_ [20] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		double sample3_ [20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };
		mcRng rng_;

	public:
		TEST_METHOD(fitAsimutDistribution)
		{
			rng_.init(33, 37);

			try
			{
				Logger::WriteMessage("In method mcSourceModelRadialPhsp::fitAsimutDistribution");

				Logger::WriteMessage("Testing sample1_");
				runSample(sample1_, true);

				Logger::WriteMessage("Testing sample2_");
				runSample(sample2_, false);

				Logger::WriteMessage("Testing sample3_");
				runSample(sample3_, false);

				Logger::WriteMessage("mcSourceModelRadialPhsp::fitAsimutDistribution tested successfully");
			}
			catch (std::out_of_range ex)
			{
				// Correct exception.
			}
			catch (const std::exception& ex)
			{
				Logger::WriteMessage((std::string("Exception in method fitAsimutDistribution: ") + ex.what()).c_str());
				Assert::Fail(L"Тест fitAsimutDistribution ", LINE_INFO());
			}
		}

	private:

		void runSample(double* sample, bool checkReproducing)
		{
			int i;
			double a, b;
			char msg[256];

			mcSourceModelRadialPhsp::fitAsimutDistribution(np_, sample, a, b);
			sprintf_s(msg, 256, "a = %f \tb = %f", a, b);
			Logger::WriteMessage(msg);

			// Vector сбора статистики симуляции апроксимированного распредедления
			vector<double> simd(20, 0);
			vector<double> pars(2);
			pars[0] = a;
			pars[1] = b;

			// Симуляция распределения
			for (i = 0; i < 100000; i++)
			{
				double f = mcSourceModelRadialPhsp::sampleAzimut(pars, rng_.rnd());
				int idx = (int)(np_ * f / PI);
				if(idx < 0 || idx >= np_)
					throw std::exception("azimut sampled angle out of range");
				simd[idx] += 1.0;
			}

			dumpNormalized(sample, np_);
			dumpNormalized(simd.data(), np_);

			if(checkReproducing)
			{
				double s1 = 0, s2 = 0;
				for (i = 0; i < np_; i++)
				{
					s1 += sample[i];
					s2 += simd[i];
				}
				s1 /= np_;
				s2 /= np_;

				double sum = 0, d2 = 0;
				for (i = 0; i < np_; i++)
				{
					double d = simd[i] / s2 - sample[i] / s1;
					d2 += d * d;
					sum += sample[i] / s1;
				}

				double dispersion = 100 * sqrt(d2) / sum;

				sprintf_s(msg, 256, "Dipersion = %0.2f percents", dispersion);
				Logger::WriteMessage(msg);

				if(dispersion > 5)
					Assert::Fail(L"Тест fitAsimutDistribution (dispersion exceed the limit)", LINE_INFO());
			}
		}

		static void dumpNormalized(const double* data, int np)
		{
			int i;
			double sum = 0;
			for (i = 0; i < np; i++)
				sum += data[i];
			sum /= np;
			if (sum == 0) sum = 1;	// вероятно нулевой вектор и будет представлен нулями

			std::strstream os;
			for (i = 0; i < np; i++)
				os << "\t" << data[i] / sum;
			os << char(0);
			//os << endl;
			Logger::WriteMessage(os.str());
		}
	};
}
