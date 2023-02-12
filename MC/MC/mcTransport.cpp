#include ".\mctransport.h"
#include "../geometry/vec3d.h"
#include "mcMedia.h"
#include "mcRng.h"
#include "mcPhysics.h"
#include "mcScore.h"
#include "mcScoreTrack.h"
#include "mcThread.h"
#include <float.h>

mcTransport::mcTransport()
	:mcObj()
	, transCutoff_phot(0)
	, transCutoff_elec(0)
	, score_(nullptr)
	, previousTransport_(nullptr)
	, nextTransport_(nullptr)
	, internalTransport_(nullptr)
	, externalTransport_(nullptr)
	, media_(nullptr)
	, defmedidx_(0)
	, defdensity_(1.0)
	, stamp_(0)
	, isMultiRegions_(false)
{
}

mcTransport::mcTransport(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x)
	:mcObj()
	, transCutoff_phot(0)
	, transCutoff_elec(0)
	, score_(nullptr)
	, previousTransport_(nullptr)
	, nextTransport_(nullptr)
	, internalTransport_(nullptr)
	, externalTransport_(nullptr)
	, media_(nullptr)
	, defmedidx_(0)
	, defdensity_(1.0)
	, isMultiRegions_(false)
{
	setPosition(orgn, z, x);
}

mcTransport::~mcTransport(void)
{
	if (score_ != nullptr) delete score_;
	for (unsigned i = 0; i < regions_.size(); i++)
		delete regions_[i];
}

#pragma region Initialization

void mcTransport::setPosition(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x)
{
	geomVector3D ax(x);
	geomVector3D az(z);
	ax.normalize();
	az.normalize();
	geomVector3D ay = az ^ ax;

	mwtot_ = geomMatrix3D::ParallelShift(-orgn.x(), -orgn.y(), -orgn.z()) *
		geomMatrix3D::BuildFromAxis(ax, ay, az);
	mttow_ = mwtot_;
	mttow_.makeInverse();
}

void mcTransport::setMediaRef(const mcMedia* media)
{
	media_ = media;
}

void mcTransport::setMediumMono(const char* mname)
{
	if (media_ == nullptr)
		throw std::exception("Media should be set before assigning transport medium");
	defmedidx_ = media_->getMediumIdx(mname);
}

void mcTransport::addRegion(mcTransport* region)
{
	// Преобразовать в мировую систему
	geomVector3D orgn = (geomVector3D(0, 0, 0) * region->MT2W()) * this->MT2W();
	geomVector3D z = (geomVector3D(0, 0, 1).transformDirection(region->MT2W())).transformDirection(this->MT2W());
	geomVector3D x = (geomVector3D(1, 0, 0).transformDirection(region->MT2W())).transformDirection(this->MT2W());
	region->setPosition(orgn, z, x);
	regions_.push_back(region);
	isMultiRegions_ = true;
}

#pragma endregion

#pragma region Transport

void mcTransport::beginTransport(mcParticle& p)
{
	if (this->media_ == nullptr)
		throw std::exception((string("Media not set in the module \"") + this->getName() + "\"").c_str());

	// Указатель на транспортный объект, в котором частица находится в данный момент
	p.transport_ = this;

	// Объекты транспорта могут ставить печать в момент											
	// BeginTransport обозначая, что частица в них находилась.
	p.regionFlags |= stamp_;

	mcParticle* particle = p.thread_->NextParticle();
	*particle = p;
	particle->p = p.p * mwtot_;
	particle->plast = p.plast * mwtot_;
	particle->u = particle->u.transformDirection(mwtot_);
	particle->dnear = 0;
	particle->mfps = HowManyMFPs(p.thread_->rng());

	particle->region.idx_ = 0;
	particle->region.medidx_ = 0;
	particle->regDensityRatio = DBL_EPSILON;

	if (score_)
		score_->ScoreFluence(*particle);

	// Транспорт в локальной системе координат
	simulate(p.thread_);
}

// GG 20140910 Если частица внутри, то не нужно переносить ее на поверхность.
// Проблема проявилась при использовании с источником C-60,
// когда частицы начинали транспортироваться внутри и код пытался перенести частицу на поверхность.
void mcTransport::beginTransportInside(mcParticle& p)
{
	if (this->media_ == nullptr)
		throw std::exception((string("Media not set in the module \"") + this->getName() + "\"").c_str());

	p.transport_ = this;
	p.regionFlags |= stamp_;

	mcParticle* particle = p.thread_->NextParticle();
	*particle = p;
	particle->p = p.p * mwtot_;
	particle->plast = p.plast * mwtot_;
	particle->u = particle->u.transformDirection(mwtot_);
	particle->dnear = 0;
	particle->mfps = HowManyMFPs(p.thread_->rng());

	particle->region.idx_ = 1;
	particle->region.medidx_ = defmedidx_;
	particle->regDensityRatio = defdensity_;

	if (score_)
		score_->ScoreFluence(*particle);

	// Транспорт в локальной системе координат
	simulate(p.thread_);
}

void mcTransport::endTransport(mcParticle* particle)
{
	// Если модуль встроет в цепочку вложенных объъектов, 
	// что определеяется наличием внутреннего или внешнего объекта или обоих,
	// то должна быть указано какая поверхность пересекается.
	// В противном случае речь о неприятной ошибки конфигурации или неправильном использования типов модулей.
	if (particle->exitSurface_ == mcParticle::temb_shit_t::Undefined &&
		(externalTransport_ != nullptr || internalTransport_ != nullptr))
		throw std::exception("All Get Distance functions of embedded transports should take care about indication wihich surface will be hitted");

	else if (particle->exitSurface_ == mcParticle::temb_shit_t::Internal && internalTransport_ == nullptr)
		throw std::exception("It should not be possible to hit internal surface if internal object does not exist or not indicated");

	mcParticle pp = **particle->thread_->CurrentParticle();
	particle->thread_->RemoveParticle();
	if (pp.ke == 0) return;

	pp.p = pp.p * mttow_;
	pp.plast = pp.plast * mttow_;
	pp.u = pp.u.transformDirection(mttow_);

	// Если пересекаем внешнюю поверхность и нет охватывающего объекта, 
	// то мы просто передаем управления стандартной однонитиевой цепочке объектов.
	if ((particle->exitSurface_ == mcParticle::temb_shit_t::External && externalTransport_ == nullptr) ||
		(externalTransport_ == nullptr && internalTransport_ == nullptr))
	{
		if (pp.u.z() < 0 && previousTransport_ != nullptr)
			previousTransport_->beginTransport(pp);
		else if (pp.u.z() >= 0 && nextTransport_ != nullptr)
			nextTransport_->beginTransport(pp);
		else if (pp.trackScore_)
			pp.trackScore_->score(particle->thread_->id(), pp.t, pp.p, pp.p + (pp.u * 100), pp.ke);
	}
	else
	{
		if (particle->exitSurface_ == mcParticle::temb_shit_t::Internal && internalTransport_ != nullptr)
			internalTransport_->beginTransportInside(pp);
		else if (particle->exitSurface_ == mcParticle::temb_shit_t::External && externalTransport_ != nullptr)
			externalTransport_->beginTransportInside(pp);
		else if (pp.trackScore_)
			pp.trackScore_->score(particle->thread_->id(), pp.t, pp.p, pp.p + (pp.u * 100), pp.ke);
	}
}

#pragma endregion

#pragma region Standard simulation

void mcTransport::simulate(mcThread* thread)
{
	mcParticle** pCurParticle = thread->CurrentParticle();
	mcParticle* particleStack = thread->ParticleStack();

	while ((*pCurParticle) >= particleStack)
	{
		mcTransport* t = (*pCurParticle)->transport_;// Указатель на транспортный объект,
		// в котором частица находится в данный момент

		// Сохраняем стартовую точку для скоринга трэка.
		// Конечная точка хранится в текущей частице.
		geomVector3D point((*pCurParticle)->p);
		enum mc_particle_t pt = (*pCurParticle)->t;
		mcRegionReference region = (*pCurParticle)->region;
		double weight = (*pCurParticle)->weight;

		if ((*pCurParticle)->mfps <= 0)
			(*pCurParticle)->mfps = HowManyMFPs(thread->rng());

		double step, edep;																					//step & edep???
		mc_move_result_t mres = t->moveParticle(*pCurParticle, step, edep);									//f12 moveParticle
		if (mres == MCMR_DISCARGE)
		{
			if (edep > 0 && t->score_ != nullptr)
				t->score_->ScorePoint(edep * weight, thread->id(), region, pt, point);
			continue;
		}

		if (edep > 0 && t->score_ != nullptr)
			t->score_->ScoreLine(edep * weight, thread->id(), region, pt, point, (*pCurParticle)->p);

		// Если частица после израсходования пути остается в транспорте, то разыгрываем взаимодействие.
		// Если частица покидает транспорт, берем следующую частицу из стэка или прерываем симуляцию.

		if (mres == MCMR_INTERUCT)
		{
			const mcPhysics* phys = t->media_->getPhysics((*pCurParticle)->t);
			const mcMedium* med = t->media_->getMedium((*pCurParticle)->t, (*pCurParticle)->region.medidx_);

			point = (*pCurParticle)->p;
			pt = (*pCurParticle)->t;
			region = (*pCurParticle)->region;
			weight = (*pCurParticle)->weight;

			// Частицы с энергией ниже критической должны быть уничтожены раньше любых расчетов транспорта.
			if (phys->Discarge(*pCurParticle, *med, edep))
			{
				if (edep > 0 && t->score_ != nullptr)
					t->score_->ScorePoint(edep * weight, thread->id(), region, pt, point);
				continue;
			}

			edep = phys->DoInterruction((*pCurParticle), med);
			if (edep > 0 && t->score_ != nullptr)
				t->score_->ScorePoint(edep * weight, thread->id(), region, pt, point);
			if ((*pCurParticle)->ke == 0)
				thread->RemoveParticle();
			else
				(*pCurParticle)->mfps = 0;

			continue;
		}
		else if (mres == MCMR_EXIT)
		{
			t->endTransport(*pCurParticle);
			continue;
		}
		else if (mres == MCMR_CONTINUE)
		{
			// По некоторым причнам частица осталась в стэке нужно продолжить ее транспорт.
			// Например, частица переместилась без изменения энергии на поверхность объекта.
			continue;
		}
		else
			throw std::exception("Unexpected partical move result in simulator");
	}
}

mc_move_result_t mcTransport::moveParticle(mcParticle* particle, double& step, double& edep)
{
	edep = 0;

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

	// Переместить частицу на поверхность, если она еще не внутри
	step = DBL_MAX;
	if (particle->region.idx_ == 0)
	{
		int iregion = 0;
		if (isMultiRegions_)
		{
			for (unsigned i = 0; i < regions_.size(); i++)
			{
				double f = regions_[i]->getDistanceOutside(*particle) + DBL_EPSILON;
				if (f < step)
				{
					iregion = i;
					step = f;
				}
			}
		}
		// Если установлен флаг типа пересекаемой поверхности как внутренней,
		// То определенно речь не о повторной возможности входа в данный объект,
		// а о преходе в следующий
		else if (particle->exitSurface_ != mcParticle::temb_shit_t::Internal)
			step = getDistanceOutside(*particle) + DBL_EPSILON;

		if (step == DBL_MAX) { // промазали, летим в следующий слой
			endTransport(particle);
			return MCMR_CONTINUE;
		}

		// При необходимости запоминаем отрезок в мировой системе координат
		if (particle->trackScore_)
			particle->trackScore_->score(particle->thread_->id(), particle->t, particle->p * mttow_, (particle->p + (particle->u * step)) * mttow_, particle->ke);

		particle->p += particle->u * step;
		particle->dnear = 0;
		particle->region.idx_ = iregion + 1;

		if (isMultiRegions_)
		{
			particle->regDensityRatio = regions_[iregion]->getDefDensity();
			particle->region.medidx_ = regions_[iregion]->getDefMedIdx();
		}
		else
		{
			particle->regDensityRatio = defdensity_;
			particle->region.medidx_ = defmedidx_;
		}

		// Возвращаемся, чтобы повторить шаг.
		// В противном случае шаг до поверхности будет включен в
		// потери энергии заряженной частицы.
		if (step > DBL_EPSILON)
			return MCMR_CONTINUE;
	}

	const mcPhysics* phys = media_->getPhysics(particle->t);
	const mcMedium* med = media_->getMedium(particle->t, particle->region.medidx_);//Параметры сред транспорта фотонов и электронов

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

	double freepath = phys->MeanFreePath(particle->ke, *med, defdensity_);
	step = freepath * particle->mfps;

	if (step < particle->dnear)
	{
		double stepRequested = step;
		edep = phys->TakeOneStep(particle, *med, step);

		// Сохраняем фрагмент трека после шага, который может меняться в процессе последнего.
		// К тому же координаты точки имеют значение после шага. Поэтому в следующих расчетах шаг отрицательный.
		if (particle->trackScore_)
			particle->trackScore_->score(particle->thread_->id(), particle->t, particle->p * mttow_, (particle->p - (particle->u * step)) * mttow_, particle->ke);

		// Шаг может быть меньше запрошенного только в случае заряженных частиц
		// и означает, что просто выполнен шаг при нерерывных потерях энергии и рассения 
		// и дискретного события не случилось.
		if (step < stepRequested)
			return MCMR_CONTINUE;
		else
			return MCMR_INTERUCT;
	}

	double dist;
	if (isMultiRegions_)
	{
		dist = regions_[particle->region.idx_ - 1]->getDistanceInside(*particle);
	}
	else
	{
		dist = getDistanceInside(*particle) + DBL_EPSILON;	// Новый трик с тем, чтобы частица чуть-чуть заступала за границу;
	}

	// HACK!!
	// По непонятным причинам частицы могут быть за пределами объекта,
	// хотя транспорт как внутри. До выяснения причин вывод сообщений о подобных событиях.
	if (dist == DBL_MAX || dist < -0.01)
	{
		cout << "Wrong particle transport inside object:" << this->getName() << endl;
		cout << "Position: " << particle->p;
		cout << "Direction: " << particle->u;
		cout << "Dnear: " << particle->dnear << endl;
		cout << "Distance = " << dist << endl;
		particle->thread_->RemoveParticle();
		return MCMR_DISCARGE;
	}

	if (step < dist)
	{
		double stepRequested = step;
		edep = phys->TakeOneStep(particle, *med, step);
		particle->mfps -= step / freepath;

		if (particle->trackScore_)
			particle->trackScore_->score(particle->thread_->id(), particle->t, particle->p * mttow_, (particle->p - (particle->u * step)) * mttow_, particle->ke);

		if (isMultiRegions_)
			particle->dnear = regions_[particle->region.idx_ - 1]->getDNearInside(particle->p);
		else
			particle->dnear = getDNearInside(particle->p);

		// !! Деликатное трудно понимаемое место.
		// Решение о событии взаимодействия принимается только если реальный шаг строго равен запрошенному шагу.
		// Следует помнить, что к данному моменту шаг определен как "step = freepath * particle->mfps;"
		// В случае заряженных частиц реальный шаг в основном определяется шагом, 
		// на котором будет потеряна определенная часть энергии (например 1%).
		// Благодаря счетчику "particle->mfps" частица будет протаскиваться через непрывное торможение
		// пока не доберется места дискретного события.

		if (step < stepRequested)
			return MCMR_CONTINUE;
		else
			return MCMR_INTERUCT;
	}
	else
	{
		// HACK! На поверхности возможно залипание, если расстояние в пределах погрешности вычислений.
		static const double epsln = DBL_EPSILON * 10;
		if (dist < epsln)
			dist = epsln;
		step = dist;
		edep = phys->TakeOneStep(particle, *med, step);

		if (particle->trackScore_)
			particle->trackScore_->score(particle->thread_->id(), particle->t, particle->p * mttow_, (particle->p - (particle->u * step)) * mttow_, particle->ke);

		particle->mfps -= step / freepath;
		if (step == dist)
		{
			// !! Для поддержки композитных модулей и модулей с невыпуклыми поверхностями 
			// (например, цилиндрические кольца) мы не устанавливаем шаг выхода а меняем индес региона на 0.
			// В этом случае данная функция будет вызвана вновь и предпринята попытка продолжить транспорт в модуле.
			// В случае отсутствия попадания частица будет, наконец, передана следующему транспорту.
			particle->region.idx_ = 0;
		}
		return MCMR_CONTINUE;
	}
}

double mcTransport::HowManyMFPs(mcRng& rng)
{
	double howMany = -log(1.0 - rng.rnd());
	return MAX(howMany, DBL_EPSILON);
}

double mcTransport::etotal() const
{
	return score_ != nullptr ? score_->etotal() : 0;
}

#pragma endregion

#pragma region Scoring

void mcTransport::setScore(mcScore* score)
{
	if (score_ != nullptr) delete score_;
	score_ = score;
	score_->setTransport(this);
}

void mcTransport::setDefaultScore(int nThreads)
{
	setScore(new mcScore(this->getName(), nThreads));
}

void mcTransport::setInternalTransport(mcTransport* t)
{
	internalTransport_ = t;

	if (t != nullptr)
	{
		// Для ускорения расчетов иницилизируем матрицы преобразования координат
		// из системы данного объекта в систему вложенного объекта.
		mtoe_ = mttow_ * t->MW2T();
	}
}

void mcTransport::setExternalTransport(mcTransport* t)
{
	externalTransport_ = t;
	t->setInternalTransport(this);
}

mcTransport* mcTransport::getInternalTransportByName(const char* name)
{
	if (strcmp(this->getName(), name) == 0)
		return this;

	auto ti = this->getInternalTransport();
	if (ti != nullptr)
	{
		// Рекурсивный поиск транспорта
		auto t = ti->getInternalTransportByName(name);
		if (t != nullptr)
			return t;
	}
	return nullptr;
}

#pragma endregion

#pragma region Геометрия (методы нужны для поддержки стандартной функции moveParticle)

double mcTransport::getDistanceInside(mcParticle& p) const
{
	implementException();
	return 0;
}

double mcTransport::getDistanceOutside(mcParticle& p) const
{
	double dist = DBL_MAX;
	if (isMultiRegions_)	// Поддержка вложений в группрвром транспорте
	{
		for (unsigned i = 0; i < regions_.size(); i++)
		{
			double f = regions_[i]->getDistanceOutside(p) + DBL_EPSILON;
			if (f < dist) dist = f;
		}
	}
	else
		implementException();
	return dist;
}

double mcTransport::getDNearInside(const geomVector3D& p) const
{
	// Разрешаем не имплементировать функцию ближайшего расстояния.
	// Возвращение нуля автоматически означает просто отказ от ускорения за счет
	// экономии на вызовах функций вычисления расстояний вдоль направления.
	return 0;
}

void mcTransport::implementException()
{
	throw std::exception("Please implement or do not use default geometry");
}

#pragma endregion

#pragma region Auxilary

void mcTransport::MoveToCoordinateSystem(const geomMatrix3D& m)
{
	mttow_ = mttow_ * m;
	mwtot_ = mttow_;
	mwtot_.makeInverse();
}

mcMediumXE* mcTransport::getParticleMedium()
{
	throw std::exception("mcTransport::getParticleMedium: method should be implemented in each transport class");
}

double mcTransport::getRegionRelDensity()
{
	throw std::exception("mcTransport::getRegionRelDensity: method should be implemented in each transport class");
}

void mcTransport::dumpVRML(ostream& os) const
{
	if (score_ != nullptr)
		score_->dumpVRML(os);
	if (isMultiRegions_)
	{
		for (unsigned i = 0; i < regions_.size(); i++)
			regions_[i]->dumpVRML(os);
	}
	else
		os << "# VRML dump for transport is not implemented" << endl;
}

void mcTransport::dump(ostream& os) const
{
	os << *(const mcObj*)this << endl;
	os << "Energy = " << etotal() << endl;
	os << "transCutoff_phot = " << transCutoff_phot << endl;
	os << "transCutoff_elec = " << transCutoff_elec << endl;
	os << "MWTOT:" << endl;
	os << mwtot_ << endl;
	os << "MTTOW:" << endl;
	os << mttow_ << endl;
	if (isMultiRegions_)
	{
		for (unsigned i = 0; i < regions_.size(); i++)
		{
			os << std::endl;
			regions_[i]->dump(os);
		}
		os << std::endl;
	}
}

#pragma endregion

#pragma region VRML utilities

void mcTransport::dumpVRMLRing(ostream& os, double r1, double r2, double z, bool normPositive, double x0, double y0) const
{
	int i, na = 24;
	double da = 6.2832 / na;
	double f1 = normPositive ? r1 : r2;
	double f2 = normPositive ? r2 : r1;

	// Преобразование из системы объекта в мировую систему 
	// с учетом смещения цилиндра в плоскости XY собственной системы координат
	geomMatrix3D mttow = geomMatrix3D::ParallelShift(x0, y0, 0) * mttow_;

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

	for (i = 0; i < na; i++) {
		geomVector3D p = geomVector3D(f1*cos(i*da), f1*sin(i*da), z) * mttow;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;
		p = geomVector3D(f2*cos(i*da), f2*sin(i*da), z) * mttow;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		if (i < na - 1) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	for (i = 0; i < na; i++) {
		os << "                " << 2 * i << ", " << 2 * i + 1 << ", " << 2 * ((i + 1) % na) + 1 << ", " << 2 * ((i + 1) % na);
		if (i < na - 1) os << ", -1,";
		os << endl;
	}

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}

void mcTransport::dumpVRMLSemiCircle(ostream& os, double r, double z, bool normPositive, const geomMatrix3D& M) const
{
	int i, na = 24;
	double da = PI / na;

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

	geomVector3D p = geomVector3D(0, 0, z) * M;
	os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;

	for (i = 0; i <= na; i++) 
	{
		double a = i*da + PI * 0.5;
		geomVector3D p = geomVector3D(r*cos(a), r*sin(a), z) * M;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		if (i < na) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	for (i = 0; i < na; i++) {
		if (normPositive)
			os << "                " << i + 1 << ", " << i + 2 << ", " << 0;
		else
			os << "                " << i + 2 << ", " << i + 1 << ", " << 0;
		if (i < na) os << ", -1,";
		os << endl;
	}

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}

void mcTransport::dumpVRMLCylinderSide(ostream& os, double r, double z1, double z2, bool normPositive, double x0, double y0) const
{
	dumpVRMLConicalCylinderSide(os, r, z1, z2, 0, normPositive, x0, y0);
}

void mcTransport::dumpVRMLConicalCylinderSide(ostream& os, double r, double z1, double z2, double f, bool normPositive, double x0, double y0) const
{
	int i, na = 24;
	double da = 6.2832 / na;
	double f1, f2, s1 = 1, s2 = 1;
	if (normPositive)
	{
		f1 = z1; f2 = z2;
		if (f != 0) s2 = (f - z2 + z1) / f;
	}
	else
	{
		f1 = z2; f2 = z1;
		if (f != 0) s1 = (f - z2 + z1) / f;
	}

	// Преобразование из системы объекта в мировую систему 
	// с учетом смещения цилиндра в плоскости XY собственной системы координат
	geomMatrix3D mttow = geomMatrix3D::ParallelShift(x0, y0, 0) * mttow_;

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

	for (i = 0; i < na; i++) {
		geomVector3D p = geomVector3D(r*s1*cos(i*da), r*s1*sin(i*da), f1) * mttow;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;
		p = geomVector3D(r*s2*cos(i*da), r*s2*sin(i*da), f2) * mttow;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		if (i < na - 1) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	for (i = 0; i < na; i++) {
		os << "                " << 2 * i << ", " << 2 * i + 1 << ", " << 2 * ((i + 1) % na) + 1 << ", " << 2 * ((i + 1) % na);
		if (i < na - 1) os << ", -1,";
		os << endl;
	}

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}

void mcTransport::dumpVRMLCylinderSemiSide(ostream& os, double r, double h, const geomMatrix3D& M) const
{
	int i, na = 24;
	double da = PI / na;

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

	for (i = 0; i <= na; i++) {
		double a = i*da + PI * 0.5;
		geomVector3D p = geomVector3D(r*cos(a), r*sin(a), h) * M;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z() << ", " << endl;
		p = geomVector3D(r*cos(a), r*sin(a), 0) * M;
		os << "                    " << p.x() << ' ' << p.y() << ' ' << p.z();
		if (i < na) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	for (i = 0; i < na; i++) {
		os << "                " << 2 * i << ", " << 2 * i + 1 << ", " << 2 * (i + 1) + 1 << ", " << 2 * (i + 1);
		if (i < na - 1) os << ", -1,";
		os << endl;
	}

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}

void mcTransport::dumpVRMLCylinder(ostream& os, double r, double z1, double z2, double x0, double y0) const
{
	dumpVRMLCylinderSide(os, r, z1, z2, false, x0, y0);
	dumpVRMLRing(os, 0, r, z1, false, x0, y0);
	dumpVRMLRing(os, 0, r, z2, true, x0, y0);
}

void mcTransport::dumpVRMLCylinderRing(ostream& os, double r1, double r2, double z1, double z2) const
{
	dumpVRMLCylinderSide(os, r1, z1, z2, true);
	dumpVRMLCylinderSide(os, r2, z1, z2, false);
	dumpVRMLRing(os, r1, r2, z1, false);
	dumpVRMLRing(os, r1, r2, z2, true);
}

void mcTransport::dumpVRMLConicalRing(ostream& os, double r1, double r2, double z1, double z2, double f) const
{
	double s = (f - z2 + z1) / f;
	dumpVRMLConicalCylinderSide(os, r1, z1, z2, f, true);
	dumpVRMLConicalCylinderSide(os, r2, z1, z2, f, false);
	dumpVRMLRing(os, r1, r2, z1, false);
	dumpVRMLRing(os, r1 * s, r2 *s, z2, true);
}

void mcTransport::dumpVRMLCylinderWithConicalHole(ostream& os, double r1, double r2, double z1, double z2, double f) const
{
	double s = (f - z2 + z1) / f;
	dumpVRMLConicalCylinderSide(os, r1, z1, z2, f, true);
	dumpVRMLCylinderSide(os, r2, z1, z2, false);
	dumpVRMLRing(os, r1, r2, z1, false);
	dumpVRMLRing(os, r1 * s, r2, z2, true);
}

void mcTransport::dumpVRMLRectangleRing(ostream& os, double x1, double x2, double y1, double y2, double d, double h) const
{
	geomVector3D p[16];
	p[0] = geomVector3D(x1 - d, y1 - d, 0) * mttow_;
	p[1] = geomVector3D(x1 - d, y2 + d, 0) * mttow_;
	p[2] = geomVector3D(x2 + d, y2 + d, 0) * mttow_;
	p[3] = geomVector3D(x2 + d, y1 - d, 0) * mttow_;
	p[4] = geomVector3D(x1, y1, 0) * mttow_;
	p[5] = geomVector3D(x1, y2, 0) * mttow_;
	p[6] = geomVector3D(x2, y2, 0) * mttow_;
	p[7] = geomVector3D(x2, y1, 0) * mttow_;
	p[8] = geomVector3D(x1 - d, y1 - d, h) * mttow_;
	p[9] = geomVector3D(x1 - d, y2 + d, h) * mttow_;
	p[10] = geomVector3D(x2 + d, y2 + d, h) * mttow_;
	p[11] = geomVector3D(x2 + d, y1 - d, h) * mttow_;
	p[12] = geomVector3D(x1, y1, h) * mttow_;
	p[13] = geomVector3D(x1, y2, h) * mttow_;
	p[14] = geomVector3D(x2, y2, h) * mttow_;
	p[15] = geomVector3D(x2, y1, h) * mttow_;

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

	for (int i = 0; i < 16; i++) {
		os << "                    " << p[i].x() << ' ' << p[i].y() << ' ' << p[i].z();
		if (i < 15) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	// Внешняя боковая поверхность
	os << "                0, 8, 9, 1, -1," << endl;
	os << "                1, 9, 10, 2, -1," << endl;
	os << "                2, 10, 11, 3, -1," << endl;
	os << "                3, 11, 8, 0, -1," << endl;

	// Внутренняя боковая поверхность
	os << "                4, 5, 13, 12, -1," << endl;
	os << "                5, 6, 14, 13, -1," << endl;
	os << "                6, 7, 15, 14, -1," << endl;
	os << "                7, 4, 12, 15, -1," << endl;

	// Нижний торец
	os << "                0, 1, 5, 4, -1," << endl;
	os << "                1, 2, 6, 5, -1," << endl;
	os << "                2, 3, 7, 6, -1," << endl;
	os << "                3, 0, 4, 7, -1," << endl;

	// Верхний торец
	os << "                8, 12, 13, 9, -1," << endl;
	os << "                9, 13, 14, 10, -1," << endl;
	os << "                10, 14, 15, 11, -1," << endl;
	os << "                11, 15, 12, 8, -1" << endl;

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}

void mcTransport::dumpVRMLPrism(ostream& os, double ax, double ay, double az) const
{
	int i = 0;
	geomVector3D p[8];

	p[i++] = geomVector3D(-0.5*ax, -0.5*ay, 0) * mttow_;
	p[i++] = geomVector3D(-0.5*ax, 0.5*ay, 0) * mttow_;
	p[i++] = geomVector3D(0.5*ax, 0.5*ay, 0) * mttow_;
	p[i++] = geomVector3D(0.5*ax, -0.5*ay, 0) * mttow_;
	p[i++] = geomVector3D(-0.5*ax, -0.5*ay, az) * mttow_;
	p[i++] = geomVector3D(-0.5*ax, 0.5*ay, az) * mttow_;
	p[i++] = geomVector3D(0.5*ax, 0.5*ay, az) * mttow_;
	p[i++] = geomVector3D(0.5*ax, -0.5*ay, az) * mttow_;

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

void mcTransport::dumpVRMLPolygonCircle(ostream& os, const std::vector<double>& pz, const std::vector<double>& pr) const
{
	if(pz.size() < 2)
		throw std::exception("dumpVRMLPolygonCircle: less thet 2 polygon points");

	for (unsigned i = 1; i < pz.size(); i++)
	{
		bool normPositive = false;
		double r1 = pr[i - 1], r2 = pr[i];
		double z1 = pz[i - 1], z2 = pz[i];
		if(z1 >= z2)
			throw std::exception("dumpVRMLPolygonCircle: Z step must be positive and not zero");

		double f = r1 == r2 ? 1e6 : r1 * (z2 - z1) / (r1 - r2);
		if (f < 0)
		{
			double r = r1; r1 = r2; r2 = r;
			double z = z1; z1 = z2; z2 = z;
			f = r1 * (z2 - z1) / (r1 - r2);
			normPositive = true;
		}
		dumpVRMLConicalCylinderSide(os, r1, z1, z2, f, normPositive);
	}

	// Торцы
	dumpVRMLRing(os, 0, pr.front(), pz.front(), false);
	dumpVRMLRing(os, 0, pr.back(), pz.back(), true);

}

#pragma endregion
