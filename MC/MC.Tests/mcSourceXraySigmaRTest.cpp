#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcSourceXraySigmaR.h"
#include "../MC/mcThread.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcSourceXraySigmaRTest)
	{
	public:

		TEST_METHOD(sample)
		{
			double sigmar = 0.1;
			double theta = 0.0;

			mcSourceXraySigmaR source("Gauss distributed electron beam", 1,
				mc_particle_t::MCP_NEGATRON, 6.0, 0.0, sigmar, theta);

			mcThread thread;
			thread.setId(0);
			mcParticle p;

			double xcount = 0;
			double x2count = 0;
			double ycount = 0;
			double y2count = 0;
			double r2count = 0;

			int i, n = 100000;
			for (i = 0; i < n; i++)
			{
				source.sample(p, &thread);
				double x = p.p.x();
				double y = p.p.y();
				xcount += x;
				x2count += x * x;
				ycount += y;
				y2count += y * y;
				r2count += x * x + y * y;
			}

			double sigmax = sqrt(x2count / n);
			double sigmay = sqrt(y2count / n);
			double sigma = sqrt(r2count / n);

			Assert::AreEqual(xcount / n, 0.0, 5.0 / sqrt(double(n)), L"sample failed", LINE_INFO());
			Assert::AreEqual(ycount / n, 0.0, 5.0 / sqrt(double(n)), L"sample failed", LINE_INFO());
			Assert::AreEqual(sigmax, sigmar * sqrt(0.5), sigmar * 5.0 / sqrt(double(n)), L"sample failed", LINE_INFO());
			Assert::AreEqual(sigmay, sigmar * sqrt(0.5), sigmar * 5.0 / sqrt(double(n)), L"sample failed", LINE_INFO());
			Assert::AreEqual(sigma, sigmar, sigmar * 5.0 / sqrt(double(n)), L"sample failed", LINE_INFO());
		}
	};
}
