#include "mcTransportSimpleSplitter.h"

mcTransportSimpleSplitter::mcTransportSimpleSplitter(
	const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, mc_particle_t ptype, int nsplit)
	:mcTransport(orgn, z, x)
	, ptype_(ptype)
	, nsplit_(nsplit)
{
}

mcTransportSimpleSplitter::~mcTransportSimpleSplitter(void)
{
}

double mcTransportSimpleSplitter::getDistanceInside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;
	return 0;
}

double mcTransportSimpleSplitter::getDistanceOutside(mcParticle& p) const
{
	p.exitSurface_ = mcParticle::temb_shit_t::External;
	return DBL_MAX;
}

void mcTransportSimpleSplitter::beginTransport(mcParticle& p)
{
	// Не трогаем стэк.
	// Частицы можно вообще транспортировать начиная с какого-то объекта.
	// Нюанс только в контроле за потоками.

	if (p.u.z() < 0 && previousTransport_ != nullptr)
		previousTransport_->beginTransport(p);
	else if (p.u.z() > 0 && nextTransport_ != nullptr)
	{
		p.weight /= nsplit_;
		for (int i = 0; i < nsplit_; i++)
			nextTransport_->beginTransport(p);
	}
}

void mcTransportSimpleSplitter::beginTransportInside(mcParticle& p)
{
	if (externalTransport_ != nullptr)
	{
		p.weight /= nsplit_;

		p.region.idx_ = 1;
		p.region.medidx_ = defmedidx_;
		p.regDensityRatio = defdensity_;

		for (int i = 0; i < nsplit_; i++)
			externalTransport_->beginTransportInside(p);
	}
}

void mcTransportSimpleSplitter::dumpVRML(ostream& os) const
{
	double a = 20;  // размера бокса
	geomVector3D p = geomVector3D(0, 0, 0.005) * mttow_;

	os << "# Slab: " << this->getName() << endl;
	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Box { size " << a << ' ' << a << ' ' << 0.01 << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
