#include "mcTransportTube3D.h"
#include "mcGeometry.h"

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

		geomVector3D X = (n.z() < n.x() && n.z() < n.y()) ? geomVector3D(0, 0, 1) :
			(n.y() < n.x()) ? geomVector3D(0, 1, 0) : geomVector3D(1, 0, 0);
		geomVector3D Y = n ^ x;
		Y.normalize();
		X = Y ^ n;
		segment.ME2SPlane = geomMatrix3D::ParallelShift(-p.x(), -p.y(), -p.z()) *
			geomMatrix3D::BuildFromAxis(X, Y, n);
		// Q: не нужно лии взять обратную матрицу ???
		
		// Аналогично система цилиндра при земене n->v.
		geomVector3D XT = (v.z() < v.x() && v.z() < v.y()) ? geomVector3D(0, 0, 1) :
			(n.y() < n.x()) ? geomVector3D(0, 1, 0) : geomVector3D(1, 0, 0);
		geomVector3D YT = v ^ x;
		YT.normalize();
		XT = YT ^ v;
		segment.ME2STube = geomMatrix3D::ParallelShift(-p.x(), -p.y(), -p.z()) *
			geomMatrix3D::BuildFromAxis(XT, YT, v);
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
			auto us = p.u * s.ME2SPlane;
			double d_side = segmentTubeDistanceInside(si, ps, us);
			double d_plane = DBL_MAX;

			// Летим назад
			if (us.z() < 0)
			{
				d_plane = -ps.z() / us.z();
				if (d_side <= d_plane) { dist += d_side; break; }
				else if (si == 0) { dist += d_plane; break; }
				else
				{
					// Переходим в предыдущий сегмент 
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
				if (d_side <= d_plane)
				{
					dist += d_side;
					break;
				}
				else if (si >= ns - 1)
				{
					dist += d_plane;
					break;
				}
				else
				{
					// Переходим в следующий сегмент 
					dist += d_plane;
					pp += p.u * d_plane;
					si++;
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

int mcTransportTube3D::segmentTubeDistanceInside(int idx, const geomVector3D& p, const geomVector3D& u) const
{
	return 0;
}

int mcTransportTube3D::segmentTubeDistanceOutside(int idx, const geomVector3D& p, const geomVector3D& u) const
{
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
