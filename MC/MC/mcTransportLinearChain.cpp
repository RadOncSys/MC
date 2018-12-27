#include "mcTransportLinearChain.h"
#include "mcGeometry.h"

mcTransportLinearChain::mcTransportLinearChain(const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x)
	: mcTransport(orgn, z, x)
{
}

mcTransportLinearChain::~mcTransportLinearChain(void)
{
	for (unsigned i = 0; i < chainTransports_.size(); i++)
		delete chainTransports_[i];
}

void mcTransportLinearChain::beginTransport(mcParticle& p)
{
	// Сюда частица попадает из одного из крайних модулей цепочки.
	// Нужно просто переправить ее внешнему модулю.
	externalTransport_->beginTransportInside(p);
}

void mcTransportLinearChain::beginTransportInside(mcParticle& p)
{
	if (p.transportNearest_ != nullptr)
	{
		auto t = p.transportNearest_;
		p.transportNearest_ = nullptr;
		p.exitSurface_ = mcParticle::temb_shit_t::External;
		t->beginTransport(p);
	}
	else
		externalTransport_->beginTransportInside(p);
}

double mcTransportLinearChain::getDistanceOutside(mcParticle& p) const
{
	// HACK! Запрещаем доступ частиц снаружи, так как по непонятным причинам возникает
	// зацикленность обращений внутрь.
	return DBL_MAX;


	p.exitSurface_ = mcParticle::temb_shit_t::External;
	double dist = DBL_MAX;
	for each (auto t in chainTransports_)
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

mcTransport* mcTransportLinearChain::getInternalTransportByName(const char* name)
{
	auto t = mcTransport::getInternalTransportByName(name);
	if (t != nullptr)
		return t;
	else if (chainTransports_.size() > 0)
	{
		for each(auto ti in chainTransports_)
		{
			auto tt = ti->getInternalTransportByName(name);
			if (tt != nullptr)
				return tt;
		}
	}
	return nullptr;
}

void mcTransportLinearChain::addTransport(mcTransport* t)
{
	chainTransports_.push_back(t);
}

void mcTransportLinearChain::completeInit()
{
	for (int i = 0; i < chainTransports_.size(); i++)
	{
		if(i > 0)
			chainTransports_[i]->setPreviousTransport(chainTransports_[i-1]);
		if(i < chainTransports_.size() - 1)
			chainTransports_[i]->setNextTransport(chainTransports_[i+1]);
	}
	chainTransports_.front()->setPreviousTransport(this);
	chainTransports_.back()->setNextTransport(this);
	
	// HACK!!! Нейтрализуем предыдущую установку.
	this->setPreviousTransport(nullptr);
	this->setNextTransport(nullptr);
}

void mcTransportLinearChain::dump(ostream& os) const
{
	mcTransport::dump(os);
	for each(auto t in chainTransports_)
		t->dump(os);
}

void mcTransportLinearChain::dumpVRML(ostream& os) const
{
	os << "# mcTransportLinearChain: " << this->getName() << endl;

	for each(auto t in chainTransports_)
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
