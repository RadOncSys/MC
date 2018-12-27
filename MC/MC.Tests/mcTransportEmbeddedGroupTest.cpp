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
#include "../MC/mcTransportEmbeddedGroup.h"
#include "../MC/mcTransportCylinder.h"
#include "../MC/mcETransportTrap.h"
#include <memory>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportEmbeddedGroupTest)
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
			std::wstring err(L"mcTransportEmbeddedGroupTest::beginTransportInside failed");

			// Сцена
			mcSourceSimpleMono source("MonoSource", 1, mc_particle_t::MCP_PHOTON, 1.0, geomVector3D(0, 0, -30.), geomVector3D(0, 0, 1.0));
			auto score = new mcScoreSphereFluence("Trap", 1);	// скоринг удаляется транспортом

			mcTransportEmbeddedGroup t_group(geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0));
			t_group.setName("Group");

			mcTransportCylinder t_cg(geomVector3D(0, 0, -50.0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), 5.0, 100);
			t_cg.setName("CG");
			t_cg.setMediaRef(&media);
			t_cg.setMediumMono("AIR700ICRU");

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
			t_group.setExternalTransport(&t_cg);
			t_group.addTransport(t_c1);
			t_group.addTransport(t_c2);

			t_group.setExternalTransport(&t_cg);
			t_cg.setExternalTransport(&t_trap);

			// Симуляции
			int i;
			source.sample(p, &thread);
			for(i=0; i<10; i++)
				t_c1->beginTransport(p);

			p.p.set(-100, 0, 0);
			for(i=0; i<10; i++)
				t_cg.beginTransport(p);

			p.p.set(0, 0, -15);
			for(i=0; i<10; i++)
				t_c1->beginTransportInside(p);

			// Поскольку транспорт реальный частицы могут достичь детектора с потерей энергии
			Assert::AreEqual(30., score->etotal(), 1.0, err.c_str(), LINE_INFO());
		}
	};
}