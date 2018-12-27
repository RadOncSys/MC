// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcMedia.h"
#include "../MC/mcThread.h"
#include "../MC/mcSourceSimpleMono.h"
#include "../MC/mcScoreSphereFluence.h"
#include "../MC/mcTransportLinearChain.h"
#include "../MC/mcTransportCylinder.h"
#include "../MC/mcETransportTrap.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportLinearChainTest)
	{
	public:

		TEST_METHOD(beginTransportInside)
		{
			mcMedia media;
			media.addName("AIR700ICRU");
			//media.initXEFromFile("../data/AcceleratorSimulator.pegs4dat");
			media.initXEFromFile("C:/Users/GennadyGorlachev/Documents/GitHub/RoissVS/MC/data/AcceleratorSimulator.pegs4dat");

			mcThread thread;
			thread.setId(0);
			mcParticle p;
			std::wstring err(L"mcTransportLinearChainTest::beginTransportInside failed");

			// Сцена
			mcSourceSimpleMono source("MonoSource", 1, mc_particle_t::MCP_PHOTON, 1.0, geomVector3D(0, 0, -60.), geomVector3D(0, 0, 1.0));
			auto score = new mcScoreSphereFluence("Trap", 1);	// скоринг удаляется транспортом

			mcTransportLinearChain t_chain(geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0));
			t_chain.setName("Chain");

			mcTransportCylinder t_cglobal(geomVector3D(0, 0, -50.0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), 15.0, 100);
			t_cglobal.setName("CylinderGlobal");
			t_cglobal.setMediaRef(&media);
			t_cglobal.setMediumMono("AIR700ICRU");

			auto t_c1 = new mcTransportCylinder(geomVector3D(0, 0, -20.0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), 5.0, 10);
			t_c1->setName("C1");
			t_c1->setMediaRef(&media);
			t_c1->setMediumMono("AIR700ICRU");

			auto t_c2 = new mcTransportCylinder(geomVector3D(0, 0, 10.0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), 5.0, 10);
			t_c2->setName("C2");
			t_c2->setMediaRef(&media);
			t_c2->setMediumMono("AIR700ICRU");

			mcETransportTrap t_trap(geomVector3D(0, 0, 50), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0));
			t_trap.setName("Trap");
			t_trap.setScore(score);

			// Иерархия объектов
			t_c1->setPreviousTransport(&t_chain);
			t_c1->setNextTransport(t_c2);
			t_c2->setPreviousTransport(t_c1);
			t_c2->setNextTransport(&t_chain);
			t_chain.setPreviousTransport(nullptr);

			t_chain.addTransport(t_c1);
			t_chain.addTransport(t_c2);

			t_cglobal.setNextTransport(&t_trap);
			t_cglobal.setPreviousTransport(&t_trap);
			t_chain.setExternalTransport(&t_cglobal);

			// Симуляции
			int i;
			source.sample(p, &thread);
			// 1
			for (i = 0; i < 10; i++)
				t_cglobal.beginTransport(p);

			// 2
			p.p.set(-30, 0, 0);
			for (i = 0; i < 10; i++)
				t_cglobal.beginTransportInside(p);

			// 3
			for (i = 0; i < 10; i++)
				t_c1->beginTransport(p);

			// 4
			p.p.set(0, 0, 0);
			for (i = 0; i < 10; i++)
				t_c2->beginTransport(p);

			// 5
			p.p.set(0, 0, 15);
			p.u.set(0, 0, -1);
			for (i = 0; i < 10; i++)
				t_c2->beginTransportInside(p);

			// 6
			p.p.set(-20, 0, -15);
			p.u.set(1, 0, 0);
			for (i = 0; i < 10; i++)
				t_cglobal.beginTransport(p);

			// 7
			p.p.set(20, 0, 15);
			p.u.set(-1, 0, 0);
			for (i = 0; i < 10; i++)
				t_cglobal.beginTransport(p);

			// Поскольку транспорт реальный частицы могут достичь детектора с потерей энергии
			Assert::AreEqual(70., score->etotal(), 1.0, err.c_str(), LINE_INFO());
		}
	};
}