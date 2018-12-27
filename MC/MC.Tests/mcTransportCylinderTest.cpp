// Radiation Oncology Monte Carlo open source project
//
// Author: [2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#include "stdafx.h"
#include "CppUnitTest.h"
#include "../MC/mcTransportCylinder.h"
#include "../geometry/vec3d.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MCTests
{
	TEST_CLASS(mcTransportCylinderTest)
	{
	public:

		TEST_METHOD(getDistanceInside)
		{
			// GG 20141220
			// Тест ситуации, реально встретившейся перед тестированием, когда было получено 
			// отрицательное расстояние заведомо большее ошибки округления.
			// Позднее выяснилось, что класс работает. Проблема внутри самого движка.
			// На самом деле координата по z была отрицательной и в этот момент выяснилось,
			// что в функции TakeOneStep для зараженных частиц не уменьшалась переменная dnear.
			// В результате заряженные частицы вываливась за пределы объекта без проверки границы.

			// Wrong particle transport inside object : Filter
			// Position : 0.313537      0.485755 - 0.0149412      1
			// Direction : 0.995603 - 0.0767558 - 0.0537046      1
			// Distance = -0.278211

			mcTransportCylinder transport(geomVector3D(0, 0, 0), geomVector3D(0, 0, 1), geomVector3D(1, 0, 0), 5.0, 0.4);

			mcParticle p;
			p.p.set(0.313537, 0.485755, 0.0149412);
			p.u.set(0.995603, -0.0767558, -0.0537046);

			double dist = transport.getDistanceInside(p);

			//Assert::IsTrue(dist > 0, L"getDistanceInside failed", LINE_INFO());
			Assert::AreEqual(0.278211, dist, 0.0001, L"getDistanceInside failed", LINE_INFO());
		}
	};
}
