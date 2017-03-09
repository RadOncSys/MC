#include "mcTransportAxialSymmetricSplitter.h"

mcTransportAxialSymmetricSplitter::mcTransportAxialSymmetricSplitter(
	const geomVector3D& orgn, const geomVector3D& z, const geomVector3D& x, mc_particle_t ptype, int nsplit)
	:mcTransport(orgn, z, x)
	, ptype_(ptype)
	, nsplit_(nsplit)
{
}

mcTransportAxialSymmetricSplitter::~mcTransportAxialSymmetricSplitter(void)
{
}

void mcTransportAxialSymmetricSplitter::beginTransport(mcParticle& p)
{
	//// Test. Отсечка заряженных частиц для проверки гипотез
	//if (p.t != mc_particle_t::MCP_PHOTON)
	//{
	//	return;
	//}

	if (p.u.z() < 0 && previousTransport_ != nullptr)
		previousTransport_->beginTransport(p);
	else if (p.u.z() > 0 && nextTransport_ != nullptr)
	{
		if (p.t == ptype_)
		{
			p.weight /= nsplit_;
			double a = 2 * PI / nsplit_;
			double sa = sin(a), ca = cos(a);
			for (int i = 0; i < nsplit_; i++)
			{
				double x = p.p.p_[0], y = p.p.p_[1];
				p.p.p_[0] = x * ca + y * sa;
				p.p.p_[1] = -x * sa + y * ca;

				x = p.u.p_[0]; y = p.u.p_[1];
				p.u.p_[0] = x * ca + y * sa;
				p.u.p_[1] = -x * sa + y * ca;

				nextTransport_->beginTransport(p);
			}
		}
		else
			nextTransport_->beginTransport(p);
	}
}

void mcTransportAxialSymmetricSplitter::dumpVRML(ostream& os) const
{
	double r = 10;  // размера бокса
	geomVector3D p = geomVector3D(0, 0, 0.005) * mttow_;

	os << "# Axial splitter: " << this->getName() << endl;
	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
	os << "  rotation 1 0 0 1.5708" << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Cylinder { " << endl;
	os << "                           radius " << r << endl;
	os << "                           height 0.01" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
