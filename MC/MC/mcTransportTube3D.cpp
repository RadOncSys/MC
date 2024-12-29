#include "mcTransportTube3D.h"
#include "mcGeometry.h"
#include "mcDefs.h"

mcTransportTube3D::mcTransportTube3D(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, 
	const std::vector<geomVector3D>& p, const std::vector<double>& r)
	: mcTransport(orgn, z, x)
{
	int ns = (int)p.size();
	segments_ = std::make_unique<std::vector<TubeSection>>(ns);

	// Настраиваем сегменты трубки.
	
	// Первый проход - устанавливаем исходные точки полигона и радиусы
	for (int i = 0; i < ns; i++)
	{
		auto& segment = segments_->at(i);
		segment.P0 = p[i];
		segment.R = r[i];
	}
	// Второй - направления осей сегментов
	for (int i = 0; i < ns - 1; i++)
	{
		auto& s0 = segments_->at(i);
		auto& s1 = segments_->at(i + 1);
		s0.V0 = s1.P0 - s0.P0;
		s0.D = s0.V0.length();
		s0.V0 /= s0.D;
	}
	segments_->at(ns - 1).V0 = segments_->at(ns - 2).V0;
	// Третий - нормали и матрицы преобразования координат
	for (int i = 0; i < ns; i++)
	{
		auto& segment = segments_->at(i);
		if (i == 0 || i == ns - 1) segment.N0 = segment.V0;
		else
			segment.N0 = (segments_->at(i - 1).V0 + segment.V0) * 0.5;

		// Система плоскости - ось Z уже определена нормалью.
		// Выбираем ось мировой системы максимально удаленную от Z 
		// и используем ее для построения ортогональной системы через векторное произведение.
		const auto& p = segment.P0;
		const auto& n = segment.N0;
		const auto& v = segment.V0;

		// Изначально транспорт разрабатывался для брахитерапии.
		// Рекомендуемая ориентация шаблона если в одной плоскости, 
		// то ось трубки (в начале) вдоль Z, а изгиб в направлении Y.
		// В таком случае при наклоне оси < 60 градусов в качестве оси X 
		// выбирается мировая ось X

		geomVector3D X = (v.z() > 0.5 || v.y() >= 0.5) ? geomVector3D(1, 0, 0) :
			(v.x() < 0.5) ? geomVector3D(0, 1, 0) : geomVector3D(0, 0, 1);

		geomVector3D Y = n ^ x;
		Y.normalize();
		X = Y ^ n;
		segment.ME2SPlane = geomMatrix3D::ParallelShift(-p.x(), -p.y(), -p.z()) *
			geomMatrix3D::BuildFromAxis(X, Y, n);
		// Q: не нужно лии взять обратную матрицу и правильная ли ориентация ???
		
		// Аналогично система цилиндра.
		Y = v ^ X;
		Y.normalize();
		X = Y ^ v;
		segment.ME2STube = geomMatrix3D::ParallelShift(-p.x(), -p.y(), -p.z()) *
			geomMatrix3D::BuildFromAxis(X, Y, v);
	}
}

double mcTransportTube3D::getDistanceInside(mcParticle& p) const
{
	int ns = (int)segments_->size();
	int si = findSegment(p.p);

	// По умолчанию 0 для борьбы с непредвиденными ошибками.
	// Это приведет просто к выходу частицы из объекта.
	double dist = 0;

	// !! Данный объект должен поддерживать вложение объектов 
	// чтобы симулровать полнотелы вариант полой трубкой.
	// Это значит, что из тела расчета расстояния 
	// до поверхности объекта нельзя делать return.
	while (true)
	{
		// На всякий случай если координаты частицы вне объекта
		// возвращаем 0 как бы обозначая что она сразу выходит из объекта.
		if (si < 0 || si >= ns - 1) break;

		// Мы потенциально перемещаем частицу по сегментам.
		// При этом используем счетчик пути.
		geomVector3D pp = p.p;

		// Расстояние до пресечения с секущей плоскостью и с боковой поверхностью (трубкой)
		double d_plane = DBL_MAX, d_side = DBL_MAX;

		while (si >= 0 && si < ns - 1)
		{
			auto& s = segments_->at(si);
			auto ps = pp * s.ME2SPlane;
			auto us = p.u.transformDirection(s.ME2SPlane);
			double d_plane = DBL_MAX;

			// Летим назад
			if (us.z() < 0)
			{
				d_plane = -ps.z() / us.z();

				auto p_tmp = pp + (p.u * d_plane);
				double r = p_tmp.lengthXY();
				if (r < s.R)
				{
					if(si == 0) { dist += d_plane; break; }
					// Переходим в предыдущий сегмент 
					dist += d_plane;
					pp += p.u * d_plane;
					si--;
				}
				else
				{
					dist += segmentTubeDistanceInside(si, ps, us);
					break;
				}
			}

			// Летим вперед
			else
			{
				s = segments_->at(si + 1);
				ps = pp * s.ME2SPlane;
				auto us = p.u.transformDirection(s.ME2SPlane);
				if (us.z() > 0)
					d_plane = -ps.z() / us.z();

				// Из-за наклона пересечения может не быть.
				// Тогда гарантировано должно быть пересечение с боковой стенкой.
				if(d_plane == DBL_MAX)
				{
					dist += segmentTubeDistanceInside(si, ps, us);
					break;
				}
				else
				{
					auto p_tmp = pp + (p.u * d_plane);
					double r = p_tmp.lengthXY();
					if (r < s.R)
					{
						if (si >= ns - 2) { dist += d_plane; break; }
						// Переходим в предыдущий сегмент 
						dist += d_plane;
						pp += p.u * d_plane;
						si--;
					}
					else
					{
						dist += segmentTubeDistanceInside(si, ps, us);
						break;
					}
				}
			}
		}
		break;
	}

	// При нулевом расстоянии до гранцы объекта есть какая-то проблема.
	// Заминаем ситуацию выводя частицу из объекта.
	if (dist != 0 && internalTransport_ != nullptr)
	{
		mcParticle pp(p);
		pp.p = pp.p * mtoe_;
		pp.u = p.u.transformDirection(mtoe_);
		double dist2 = internalTransport_->getDistanceOutside(pp);
		if (dist2 < dist)
		{
			p.transportNearest_ = pp.transportNearest_;
			p.exitSurface_ = mcParticle::temb_shit_t::Internal;
			return dist2;
		}
		else
		{
			p.exitSurface_ = mcParticle::temb_shit_t::External;
			return dist;
		}
	}
	else
	{
		p.exitSurface_ = mcParticle::temb_shit_t::External;
		return dist;
	}
}

double mcTransportTube3D::getDistanceOutside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;

	int ns = (int)segments_->size();
	int si = findSegment(p.p);
	double dist = 0;

	// Мы потенциально перемещаем частицу по сегментам.
	geomVector3D pp = p.p;

	// Расстояние до пресечения с секущей плоскостью и с боковой поверхностью (трубкой)
	double d_plane = DBL_MAX, d_side = DBL_MAX;

	while (true)
	{
		auto& s = segments_->at(si);
		auto us = p.u * s.ME2SPlane;

		// Заведомо летим от при том что столкновений так и не обнаружено.
		if ((si < 0 && us.z() <= 0) || (si >= ns - 1 && us.z() >= 0)) 
			return DBL_MAX;

		auto ps = pp * s.ME2SPlane;
		double d_side = segmentTubeDistanceOutside(si, ps, us);
		double d_plane = DBL_MAX;

		// Летим назад
		if (us.z() < 0)
		{
			d_plane = -ps.z() / us.z();
			if (d_side <= d_plane) 
				return dist + d_side;
			else
			{
				dist += d_plane;
				pp += p.u * d_plane;
				si--;
			}
		}

		// Летим вперед
		else
		{
			s = segments_->at(si + 1);
			ps = pp * s.ME2SPlane;
			us = p.u * s.ME2SPlane;

			if (us.z() > 0)
				d_plane = -ps.z() / us.z();
			// Не столкнулись ни с одной плоскостью сегмента, значит летим в никуда,
			// если только не столкнулись с трубкой в пределах сегмента
			else
			{
				if (d_side == DBL_MAX || si >= ns - 1) return DBL_MAX;
				else
				{
					auto ppss = (pp + p.u * d_side) * segments_->at(si + 1).ME2SPlane;
					if (ppss.z() < 0) return dist + d_side;
					else return DBL_MAX;
				}
			}

			if (d_side <= d_plane)
				return dist + d_side;
			else
			{
				dist += d_plane;
				pp += p.u * d_plane;
				si++;
			}

			// Если сегмент последний то столкновений уже не будет
			if (si >= ns - 1)
				return DBL_MAX;
		}
	}

	return DBL_MAX;
}

int mcTransportTube3D::findSegment(const geomVector3D& p) const
{
	for (int i = 0; i < segments_->size(); i++)
	{
		geomVector3D pplane = p * segments_->at(i).ME2SPlane;
		if (pplane.z() < 0 && i == 0) return -1;
		else if (pplane.z() >= 0) return i;
	}
	return (int)segments_->size();
}

double mcTransportTube3D::segmentTubeDistanceInside(int idx, const geomVector3D& p, const geomVector3D& u) const
{
	// Попытки решить аналитически ни к чему не привели.
	// Поэтому итерационное решение.
	// В основе поиск угла поворота плоскости вокруг оси Z сегмента до нахождения положения, 
	// когда в одной точке пересекаются траектория плоскость и боковая поверхность сегмента.

	// Игнорируем нехорошие ситуции вырожденности, так как по условию
	// использования пересечение с боковой стенкой гарантировано.

	// Расчеты проводятся в системе координат ME2STube сегмента.
	// Угол t отсчитывается от оси X в сторону оси Y.

	auto& s = segments_->at(idx);
	auto& s1 = segments_->at(idx + 1);

	// Траектория частицы в системе трубки
	auto pt = p * s.ME2STube;
	auto ut = u.transformDirection(s.ME2STube);

	// Направо или налево летит частица если смотреть на нее из центра координат?
	bool isRight = (pt ^ ut).z() > 0;

	double x = ut.x(), y = ut.y();
	if (abs(ut.z()) >= 1.0 - MINDELTA) { x = pt.x(), y = pt.y(); }

	double f = sqrt(x * x + y * y);
	double sint = -y / f;
	double cost = x / f;
	if (!isRight) { sint = -sint; cost = -cost; }

	double dist = 0, dist_prev = 0;
	double rp_prev =0;
	double rs_prev =0;
	double sinda = 0.1; // стартовый шаг по углу в радианах (порядка 5 градусов)

	// Итеративный поиск ограниченный 20-ю шагами
	int count = 0;
	for (; count < 20; count++)
	{
		// Нормаль к секущей плоскости
		auto NP = geomVector3D(sint, -cost, 0) ^ geomVector3D(0, 0, 1);

		// Точка пересечения траектории с секущей плоскостью
		dist = -(pt * NP) / (ut * NP);
		auto pc = pt + (ut * dist);

		// Две точки пересечения торцевых кругов с найденной плоскостью. 
		geomVector3D c0 = (pc ^ s.N0) ^ s.N0;
		c0.normalize();
		c0 = c0 * s1.R;
		geomVector3D c1 = ((pc - geomVector3D(0, 0, s.D)) ^ s1.N0) ^ s1.N0;
		c1.normalize();
		c1 = geomVector3D(0, 0, s.D) + (c1 * s1.R);

		// Вычисляем точку поверхности при координате pc.Z в секущей плоскости.
		auto ps = c0 + ((c1 - c0) * ((pc.z() - c0.z()) / (c1.z() - c0.z())));

		// Расстояния до оси
		double rp = pc.lengthXY();	// particle
		double rs = ps.lengthXY();	// surface

		// Первый шаг отличается от остальных тем, что только определяется начальное положение
		// секущей плоскости и угол отклонения может достигать 90 градусов.
		if (count == 0)
		{
			// Пересечение траектории с секущей плоскостью внутри объекта
			if (rp < rs)
			{
				// В приближении круга прогнозируем пересечение траектории с поверхностью.
				// Новая секущая плоскость проходит через спрогнозированную точку,
				// отклоняющуюся от текущей на угол da
				double cos_da = rp / rs;
				double sin_da = sqrt(1 - cos_da * cos_da);
				if (!isRight) sin_da = -sin_da;
				// Поворот новой секущей плоскости
				double a = sint;
				sint = a * cos_da + cost * sin_da;
				cost = cost * cos_da - a * sin_da;
			}
			// Пересечение за пределами объекта
			else
			{
				// Персечение где-то между положением частицы и точкой пересечения траектории с плоскостью.
				// Строим биссектрису и снова определяем в какой половине пересечение.
				// Оцениваем эквивалентный радиус как rs = rp + (rp - rs) и далае аналогично предыдущему.
				double cos_da = rs / rp;
				double sin_da = sqrt(1 - cos_da * cos_da);
				if (!isRight) sin_da = -sin_da;
				double a = sint;
				sint = a * cos_da + cost * sin_da;
				cost = cost * cos_da - a * sin_da;
			}
		}
		// На втором шаге мы только начинаем прощупывать окрестности
		else if(count == 1)
		{
			// Поворачиваем на da и на следующем шаге увидим как изменились радиусы
			double cos_da = sqrt(1 - sinda * sinda);
			double a = sint;
			sint = a * cos_da + cost * sinda;
			cost = cost * cos_da - a * sinda;
		}
		else
		{
			// Линейно интерполируем новый поворот чтобы получить нулевую разницу между радиусами
			double dr = rs - rp;
			double dr_prev = rs_prev - rp_prev;
			if (dr_prev != dr)
				sinda *= -(1 + dr_prev / (dr - dr_prev));
			double cos_da = sqrt(1 - sinda * sinda);
			double a = sint;
			sint = a * cos_da + cost * sinda;
			cost = cost * cos_da - a * sinda;
		}

		// Интерполяция поворота плоскости начинается со второго круга, 
		// но только на следующем будет получен результат.
		if (count > 2 && abs(dist_prev - dist) < MINDELTA)
			break;

		rp_prev = rp;
		rs_prev = rs;
		dist_prev = dist;
	}

	// Если не нашли хорошее решение возвращаем 0 чтобы частица просто покинула объекта
	if (count == 20)
		return 0;
	else
		return dist;
}

double mcTransportTube3D::segmentTubeDistanceOutside(int idx, const geomVector3D& p, const geomVector3D& u) const
{
	auto& s = segments_->at(idx);
	auto& s1 = segments_->at(idx + 1);

	// Траектория частицы в системе трубки
	auto pt = p * s.ME2STube;
	auto ut = u.transformDirection(s.ME2STube);

	// Направо или налево летит частица если смотреть на нее из центра координат?
	bool isRight = (pt ^ ut).z() > 0;

	double x = ut.x(), y = ut.y();
	if (abs(ut.z()) >= 1.0 - MINDELTA) { x = pt.x(), y = pt.y(); }

	double f = sqrt(x * x + y * y);
	double sint = -y / f;
	double cost = x / f;
	if (!isRight) { sint = -sint; cost = -cost; }

	double dist = 0, dist_prev = 0;
	double rp_prev = 0;
	double rs_prev = 0;
	double sinda = 0.1; // стартовый шаг по углу в радианах (порядка 5 градусов)

	// Итеративный поиск ограниченный 20-ю шагами
	int count = 0;
	for (; count < 20; count++)
	{
		// Нормаль к секущей плоскости
		auto NP = geomVector3D(sint, -cost, 0) ^ geomVector3D(0, 0, 1);

		// Точка пересечения траектории с секущей плоскостью
		dist = -(pt * NP) / (ut * NP);
		auto pc = pt + (ut * dist);

		// Две точки пересечения торцевых кругов с найденной плоскостью. 
		geomVector3D c0 = (pc ^ s.N0) ^ s.N0;
		c0.normalize();
		c0 = c0 * s1.R;
		geomVector3D c1 = ((pc - geomVector3D(0, 0, s.D)) ^ s1.N0) ^ s1.N0;
		c1.normalize();
		c1 = geomVector3D(0, 0, s.D) + (c1 * s1.R);

		// Вычисляем точку поверхности при координате pc.Z в секущей плоскости.
		auto ps = c0 + ((c1 - c0) * ((pc.z() - c0.z()) / (c1.z() - c0.z())));

		// Расстояния до оси
		double rp = pc.lengthXY();	// particle
		double rs = ps.lengthXY();	// surface

		// Первый шаг отличается от остальных тем, что только определяется начальное положение
		// секущей плоскости и угол отклонения может достигать 90 градусов.
		if (count == 0)
		{
			// Пересечение траектории с секущей плоскостью внутри объекта
			if (rp < rs)
			{
				// В приближении круга прогнозируем пересечение траектории с поверхностью.
				// Новая секущая плоскость проходит через спрогнозированную точку,
				// отклоняющуюся от текущей на угол da
				double cos_da = rp / rs;
				double sin_da = sqrt(1 - cos_da * cos_da);
				if (!isRight) sin_da = -sin_da;
				// Поворот новой секущей плоскости
				double a = sint;
				sint = a * cos_da + cost * sin_da;
				cost = cost * cos_da - a * sin_da;
			}
			// Пересечение за пределами объекта
			else
			{
				// Персечение где-то между положением частицы и точкой пересечения траектории с плоскостью.
				// Строим биссектрису и снова определяем в какой половине пересечение.
				// Оцениваем эквивалентный радиус как rs = rp + (rp - rs) и далае аналогично предыдущему.
				double cos_da = rs / rp;
				double sin_da = sqrt(1 - cos_da * cos_da);
				if (!isRight) sin_da = -sin_da;
				double a = sint;
				sint = a * cos_da + cost * sin_da;
				cost = cost * cos_da - a * sin_da;
			}
		}
		// На втором шаге мы только начинаем прощупывать окрестности
		else if (count == 1)
		{
			// Поворачиваем на da и на следующем шаге увидим как изменились радиусы
			double cos_da = sqrt(1 - sinda * sinda);
			double a = sint;
			sint = a * cos_da + cost * sinda;
			cost = cost * cos_da - a * sinda;
		}
		else
		{
			// Линейно интерполируем новый поворот чтобы получить нулевую разницу между радиусами
			double dr = rs - rp;
			double dr_prev = rs_prev - rp_prev;
			if (dr_prev != dr)
				sinda *= -(1 + dr_prev / (dr - dr_prev));
			double cos_da = sqrt(1 - sinda * sinda);
			double a = sint;
			sint = a * cos_da + cost * sinda;
			cost = cost * cos_da - a * sinda;
		}

		// Интерполяция поворота плоскости начинается со второго круга, 
		// но только на следующем будет получен результат.
		if (count > 2 && abs(dist_prev - dist) < MINDELTA)
			break;

		rp_prev = rp;
		rs_prev = rs;
		dist_prev = dist;
	}

	// Если не нашли хорошее решение возвращаем 0 чтобы частица просто покинула объекта
	if (count == 20)
		return 0;
	else
		return dist;





	return 0;
}

void mcTransportTube3D::dump(ostream& os) const
{
	mcTransport::dump(os);
	os << "Polygon:\t" << endl;
	unsigned i;
	//for (i = 0; i < pz_.size(); i++) os << pz_[i] << "\t";
	//os << endl;
	//for (i = 0; i < pr_.size(); i++) os << pr_[i] << "\t";
	os << endl << endl;
}

void mcTransportTube3D::dumpVRML(ostream& os) const
{
	os << "# mcTransportTube3D: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	//dumpVRMLPolygonCircle(os, pz_, pr_);

	os << "  ]" << endl;
	os << "}" << endl;
}
