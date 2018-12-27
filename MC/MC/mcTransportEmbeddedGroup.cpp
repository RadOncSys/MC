#include "mcTransportEmbeddedGroup.h"
#include "mcGeometry.h"

mcTransportEmbeddedGroup::mcTransportEmbeddedGroup(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x)
	: mcTransport(orgn, z, x)
{
}

mcTransportEmbeddedGroup::~mcTransportEmbeddedGroup(void)
{
	for (unsigned i = 0; i < embeddedTransports_.size(); i++)
		delete embeddedTransports_[i];
}

void mcTransportEmbeddedGroup::beginTransportInside(mcParticle& p)
{
	// Кому передать управление зависит от того кто послал частицу.
	// Если это был родитель, то частица должна быть передана одному из объектов группы.
	// Чтобы не повторять вычисление расстояния до объекта в частице предусмотрен
	// указатель еще на один транспорт, который устанавливается при предыдущем измерении

	if (p.transportNearest_ != nullptr)
	{
		auto t = p.transportNearest_;
		p.transportNearest_ = nullptr;
		p.exitSurface_ = mcParticle::temb_shit_t::External;
		t->beginTransport(p);
	}
	else
	{
		externalTransport_->beginTransportInside(p);
	}
}

void mcTransportEmbeddedGroup::endTransport(mcParticle* particle)
{
	__super::endTransport(particle);
}

void mcTransportEmbeddedGroup::addTransport(mcTransport* t)
{
	t->setExternalTransport(this);
	embeddedTransports_.push_back(t);

	// HACK!!! Нейтрализуем предыдущую установку.
	this->setInternalTransport(nullptr);
}

double mcTransportEmbeddedGroup::getDistanceInside(mcParticle& p) const
{
	return 0;
}

double mcTransportEmbeddedGroup::getDistanceOutside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;
	double dist = DBL_MAX;
	for each (auto t in embeddedTransports_)
	{
		// Нужно переводить в систему объекта
		geomVector3D p_old = p.p;
		geomVector3D u_old = p.u;
		p.p = p.p * t->MW2T();
		p.u = p.u.transformDirection(t->MW2T());

		double d = t->getDistanceOutside(p);
		if (d < dist)
		{
			dist = d;
			p.transportNearest_ = t;
		}

		p.p = p_old;
		p.u = u_old;
	}
	return dist;
}

mcTransport* mcTransportEmbeddedGroup::getInternalTransportByName(const char* name)
{
	auto t = mcTransport::getInternalTransportByName(name);
	if (t != nullptr)
		return t;
	else if (embeddedTransports_.size() > 0)
	{
		for each(auto ti in embeddedTransports_)
		{
			auto tt = ti->getInternalTransportByName(name);
			if (tt != nullptr)
				return tt;
		}
	}
	return nullptr;
}

void mcTransportEmbeddedGroup::dump(ostream& os) const
{
	mcTransport::dump(os);
	for each(auto t in embeddedTransports_)
		t->dump(os);
}

void mcTransportEmbeddedGroup::dumpVRML(ostream& os) const
{
	os << "# mcTransportEmbeddedGroup: " << this->getName() << endl;

	for each(auto t in embeddedTransports_)
	{
		t->dumpVRML(os);
		mcScore* s = t->getScore();
		if (s) s->dumpVRML(os);

		// Вложения
		mcTransport* ti = t;
		if (ti != nullptr) ti = ti->getInternalTransport();
		while (ti)
		{
			ti->dumpVRML(os);
			mcScore* ss = ti->getScore();
			if (s) ss->dumpVRML(os);
			ti = ti->getInternalTransport();
		}
	}
}
