#include "mcTransportGridFilter.h"
#include "mcGeometry.h"
#include "mcMedia.h"
#include "mcPhysics.h"
#include "mcThread.h"
#include <float.h>

#define GEOM_EPSILON	1E-6

mcTransportGridFilter::mcTransportGridFilter(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
	int nx, int ny, int nz, double psx, double psy, double psz)
	: mcTransportPrism(orgn, vz, vx, psx * nx, psy * ny, psz * nz)
	, nx_(nx), ny_(ny), nz_(nz), psx_(psx), psy_(psy), psz_(psz)
{
	x0_ = -0.5 * psx_ * nx_;
	y0_ = -0.5 * psy_ * ny_;
	z0_ = 0;

	bx_.resize(nz_, 0);
	by_.resize(nz_, 0);
}

mcTransportGridFilter::~mcTransportGridFilter(void)
{
}

int mcTransportGridFilter::getIdxAtPoint(const geomVector3D& p, short* pgidx, bool& isInBrick) const
{
	double dx = p.x() - x0_, dy = p.y() - y0_, dz = p.z() - z0_;
	if (dx < 0 || dy < 0 || dz < 0 || dx >= ax_ || dy >= ay_ || dz >= az_)
		return -1;
	else
	{
		int i = int(dx / psx_), j = int(dy / psy_), k = int(dz / psz_);
		pgidx[0] = i; pgidx[1] = j; pgidx[2] = k;

		double ddx = dx - (i + 0.5) * psx_;
		double ddy = dy - (j + 0.5) * psy_;

		if (abs(ddx) < bx_[k] * 0.5 && abs(ddy) < by_[k] * 0.5)
			isInBrick = true;
		else
			isInBrick = false;

		return (k * ny_ + j) * nx_ + i;
	}
}

double mcTransportGridFilter::getDistanceInsideVoxel(const mcParticle& particle, 
	short* gidxNext, int& idx, bool& isHitCell)
{
	int k = particle.region.gidx_[2];
	double bx = 0.5 * bx_[k], by = 0.5 * by_[k];
	idx = particle.region.idx_;
	gidxNext[0] = particle.region.gidx_[0];
	gidxNext[1] = particle.region.gidx_[1];
	gidxNext[2] = k;

	double dx = DBL_MAX, dy = DBL_MAX, dz = DBL_MAX;
	double vx = particle.u.x(), vy = particle.u.y(), vz = particle.u.z();

	// Частица внутри брикета
	if (particle.region.subidx_ > 0)
	{
		// Координаты частицы относительно границы брикета
		double x = particle.p.x() - ((x0_ + psx_ * (particle.region.gidx_[0] + 0.5)) - bx);
		double y = particle.p.y() - ((y0_ + psy_ * (particle.region.gidx_[1] + 0.5)) - by);
		double z = particle.p.z() - (z0_ + psz_ * k);

		if (vx < 0) dx = -x / vx;
		else if (vx > 0) dx = (2 * bx - x) / vx;

		if (vy < 0) dy = -y / vy;
		else if (vy > 0) dy = (2 * by - y) / vy;

		if (vz < 0) dz = -z / vz;
		else if (vz > 0) dz = (psz_ - z) / vz;

		if (dx < dy && dx < dz)
		{
			gidxNext[0] = particle.region.gidx_[0] + (vx < 0 ? -1 : 1);
			if (gidxNext[0] < 0 || gidxNext[0] >= (short)nx_) idx = -1;
			else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
			isHitCell = false;
			return dx;
		}
		else if (dy < dx && dy < dz)
		{
			gidxNext[1] = particle.region.gidx_[1] + (vy < 0 ? -1 : 1);
			if (gidxNext[1] < 0 || gidxNext[1] >= (short)ny_) idx = -1;
			else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
			isHitCell = false;
			return dy;
		}
		else
		{
			gidxNext[2] = k + (vz < 0 ? -1 : 1);
			if (gidxNext[2] < 0 || gidxNext[2] >= (short)nz_) idx = -1;
			else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
			isHitCell = true;
			return dz;
		}
	}

	// Частица за пределами брикета
	else
	{
		// Пересечение с ячейкой
		// Координаты частицы относительно границы ячейки
		double x = particle.p.x() - (x0_ + psx_ * particle.region.gidx_[0]);
		double y = particle.p.y() - (y0_ + psy_ * particle.region.gidx_[1]);
		double z = particle.p.z() - (z0_ + psz_ * k);

		if (vx < 0) dx = -x / vx;
		else if (vx > 0) dx = (psx_ - x) / vx;

		if (vy < 0) dy = -y / vy;
		else if (vy > 0) dy = (psy_ - y) / vy;

		if (vz < 0) dz = -z / vz;
		else if (vz > 0) dz = (psz_ - z) / vz;

		// Пересечения с брикетом.
		double db = mcGeometry::getDistanceToRectanglePipeOutside(
			geomVector3D(x, y, z), particle.u, psx_ / 2 - bx, psx_ / 2 + bx, psy_ / 2 - by, psy_ / 2 + by);

		if (db < dx && db < dy && db < dz)
		{
			isHitCell = false;
			return db;
		}

		else if (dx < dy && dx < dz)
		{
			gidxNext[0] = particle.region.gidx_[0] + (vx < 0 ? -1 : 1);
			if (gidxNext[0] < 0 || gidxNext[0] >= (short)nx_) idx = -1;
			else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
			isHitCell = true;
			return dx;
		}
		else if (dy < dx && dy < dz)
		{
			gidxNext[1] = particle.region.gidx_[1] + (vy < 0 ? -1 : 1);
			if (gidxNext[1] < 0 || gidxNext[1] >= (short)ny_) idx = -1;
			else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
			isHitCell = true;
			return dy;
		}
		else
		{
			gidxNext[2] = k + (vz < 0 ? -1 : 1);
			if (gidxNext[2] < 0 || gidxNext[2] >= (short)nz_) idx = -1;
			else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
			isHitCell = true;
			return dz;
		}
	}
}

void mcTransportGridFilter::beginTransport(mcParticle& p)
{
	// Указатель на транспортный объект, в котором частица находится в данный момент
	p.transport_ = this;

	mcParticle* particle = p.thread_->NextParticle();
	*particle = p;
	particle->p = p.p * mwtot_;
	particle->plast = p.plast * mwtot_;
	particle->u = particle->u.transformDirection(mwtot_);
	particle->dnear = 0;
	particle->mfps = HowManyMFPs(p.thread_->rng());

	// Переместить частицу на поверхность, если она еще не внутри
	bool isInBrick;
	int idx = getIdxAtPoint(particle->p, particle->region.gidx_, isInBrick);
	if (idx < 0)
	{
		double f = getDistanceOutside(*particle);
		if (f == DBL_MAX)
			return endTransport(particle);
		particle->p += particle->u * (f + FLT_EPSILON);
		idx = getIdxAtPoint(particle->p, particle->region.gidx_, isInBrick);
	}
	particle->region.idx_ = idx;
	particle->region.medidx_ = isInBrick ? inBrickIdx_ : outBrickIdx_;
	particle->region.subidx_ = isInBrick ? 1 : 0;
	particle->regDensityRatio = 1.0;

	// Транспорт в локальной системе координат
	simulate(p.thread_);
}

void mcTransportGridFilter::beginTransportInside(mcParticle& p)
{
	beginTransport(p);
}

mc_move_result_t mcTransportGridFilter::moveParticle(mcParticle* particle, double& step, double& edep)
{
	edep = 0;
	step = 0;

	// Важный момент! 
	// Если частица находится за пределами (на что указывает индекс региона), то возвращаемся с метокй о выходе.
	if (particle->region.idx_ < 0)
	{
		particle->exitSurface_ = mcParticle::temb_shit_t::External;
		return MCMR_EXIT;
	}

	// HACK!!
	// По непонятным причинам координаты частицы могут быть абсурдными.
	// Удалаяем такие частицы
	if (_isnan(particle->p.x()) != 0)
	{
		//cout << "Non number position or direction in object: " << this->getName() << endl;
		cout << "Non number position in object: " << this->getName() << endl;
		cout << "Position: " << particle->p;
		cout << "Direction: " << particle->u;
		particle->thread_->RemoveParticle();
		return MCMR_DISCARGE;
	}

	const mcPhysics* phys = media_->getPhysics(particle->t);

	//Параметры сред транспорта фотонов и электронов
	const mcMedium* med = media_->getMedium(particle->t, particle->region.medidx_);

	// Частицы с энергией ниже критической должны быть уничтожены раньше любых расчетов транспорта.
	if (phys->Discarge(particle, *med, edep))
		return MCMR_DISCARGE;

	// Hack!!! GG 20171030
	if (_isnan(particle->ke) != 0)
	{
		cout << "Non number energy: " << this->getName() << endl;
		cout << "Position: " << particle->p;
		cout << "Direction: " << particle->u;
		particle->thread_->RemoveParticle();
		return MCMR_DISCARGE;
	}

	double freepath = phys->MeanFreePath(particle->ke, *med, particle->regDensityRatio);
	step = freepath * particle->mfps;

	// Расстояние до границы ячейки в направлении частицы.
	// Если потом мы сделаем шаг до границы, то должны переместить частиу в следующую ячейку.
	// Именно в это момент мы пометим частицу как находящуюся на грани
	short gidxNext[3];
	int idx;
	bool isHitCell;	// true - будет пересечена граница ячейки, false - граница брикета.
	double dist = getDistanceInsideVoxel(*particle, gidxNext, idx, isHitCell) + GEOM_EPSILON;

	// Игнорируем ускорение за счет оценки расстояния до границы
	particle->dnear = 0;

	if (step < dist)
	{
		double stepRequested = step;
		edep = phys->TakeOneStep(particle, *med, step);

		if (step < stepRequested)
			return MCMR_CONTINUE;
		else
			return MCMR_INTERUCT;
	}
	else
	{
		// HACK! На поверхности возможно залипание, если расстояние в пределах погрешности вычислений.
		if (dist < GEOM_EPSILON)
			dist = GEOM_EPSILON;
		step = dist;
		edep = phys->TakeOneStep(particle, *med, step);

		particle->mfps -= step / freepath;
		if (step == dist)
		{
			// Добрались до границы воксела. Нужно перенастороить частицу на новый.
			particle->region.idx_ = idx;
			if (idx >= 0)
			{
				// Если пересекаем границу ячейки, то уходим в другую и 
				// нужно определять где именно частица окажется.
				// проще всего это сделать с помощью уже имеющейся 
				// функции определения где находится частица по ее координатам.
				if (isHitCell)
				{
					bool isInBrick;
					idx = getIdxAtPoint(particle->p, particle->region.gidx_, isInBrick);
					particle->region.idx_ = idx;
					particle->region.medidx_ = isInBrick ? inBrickIdx_ : outBrickIdx_;
					particle->region.subidx_ = isInBrick ? 1 : 0;
				}
				else
				{
					if (particle->region.subidx_ == 1)
					{
						particle->region.medidx_ = outBrickIdx_;
						particle->region.subidx_ = 0;
					}
					else
					{
						particle->region.medidx_ = inBrickIdx_;
						particle->region.subidx_ = 1;
					}
				}
			}
		}
		return MCMR_CONTINUE;
	}
}

void mcTransportGridFilter::dump(ostream& os) const
{
	__super::dump(os);
}

void mcTransportGridFilter::dumpVRML(ostream& os) const
{
	double az = psz_;
	for (int k = 0; k < nz_; k++)
	{
		double z0 = z0_ + psz_ * k;
		double ax = bx_[k];
		double ay = by_[k];
		for (int j = 0; j < ny_;j++)
		{
			double y0 = y0_ + psy_ * (j + 0.5) - 0.5 * ay;
			for (int ii = 0; ii < nx_; ii++)
			{
				double x0 = x0_ + psx_ * (ii + 0.5) - 0.5 * ax;

				// Каждый брикет изображаем как отдельный объект
				geomVector3D p[8];

				int i = 0;
				p[i++] = geomVector3D(x0, y0, z0) * mttow_;
				p[i++] = geomVector3D(x0, y0 + ay, z0) * mttow_;
				p[i++] = geomVector3D(x0 + ax, y0 + ay, z0) * mttow_;
				p[i++] = geomVector3D(x0 + ax, y0, z0) * mttow_;
				p[i++] = geomVector3D(x0, y0, z0 + az) * mttow_;
				p[i++] = geomVector3D(x0, y0 + ay, z0 + az) * mttow_;
				p[i++] = geomVector3D(x0 + ax, y0 + ay, z0 + az) * mttow_;
				p[i++] = geomVector3D(x0 + ax, y0, z0 + az) * mttow_;

				os << "    Transform {" << endl;
				os << "      children Shape {" << endl;
				os << "        appearance Appearance {" << endl;
				os << "          material Material {" << endl;
				os << "            diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
				os << "            transparency " << transparancy_ << endl;
				os << "          }" << endl;
				os << "        }" << endl;
				os << "        geometry IndexedFaceSet {" << endl;
				os << "            coord Coordinate {" << endl;
				os << "                point [" << endl;

				for (i = 0; i < 8; i++) {
					os << "                    " << p[i].x() << ' ' << p[i].y() << ' ' << p[i].z();
					if (i < 7) os << ", ";
					os << endl;
				}

				os << "                ]" << endl;
				os << "            }" << endl;
				os << "            coordIndex [" << endl;

				os << "                0, 1, 2, 3, -1," << endl;
				os << "                0, 4, 5, 1, -1," << endl;
				os << "                1, 5, 6, 2, -1," << endl;
				os << "                2, 6, 7, 3, -1," << endl;
				os << "                3, 7, 4, 0, -1," << endl;
				os << "                4, 7, 6, 5, -1" << endl;

				os << "            ]" << endl;
				os << "        }" << endl;
				os << "      }" << endl;
				os << "    }" << endl;
			}
		}
	}
}
