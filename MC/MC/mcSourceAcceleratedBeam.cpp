#include "mcSourceAcceleratedBeam.h"
#include "mcDefs.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "../geometry/text.h"

mcSourceAcceleratedBeam::mcSourceAcceleratedBeam(const char* name, int nThreads, double z)
	:mcSource(name, nThreads)
	, z_(z)
	, nparticles_(0)
	, idx_(nThreads, 0)
{
}

void mcSourceAcceleratedBeam::sample(mcParticle& p, mcThread* thread)
{
	int i = (idx_[thread->id()]++) % nparticles_;
	p = particles_[i];
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke * p.weight;
}

void mcSourceAcceleratedBeam::loadData(istream & is)
{
	mcParticle particle;
	particle.t = mc_particle_t::MCP_NEGATRON;
	particle.q = -1;
	particle.thread_ = nullptr;
	particle.trackScore_ = nullptr;

	vector<double> data(8,0);
	string line;
	std::getline(is, line, '\n');

	while (!is.fail())
	{
		std::getline(is, line, '\n');
		if (is.fail() || line.length() < 30)
			break;
		if (GetFloatArray(line, data) != 8)
			throw std::exception("mcSourceAcceleratedBeam::loadData: 8 parameters expected");

		particle.ke = data[6];
		particle.p.set(data[0], data[1], data[2] + z_);
		particle.u.set(data[3], data[4], data[5]);
		particle.weight = data[7];

		particles_.push_back(particle);
	}
	nparticles_ = (int)particles_.size();
}

void mcSourceAcceleratedBeam::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	double r = 0.1, rr = 0.4, h = 1.0, L = 2.0;

	os << "# Source: " << name_ << endl;

	geomVector3D p(0, 0, z_ - L*0.5 - h);

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
	os << "                           height " << L << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;

	p(2) = z_ - h*0.5;

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
	os << "      geometry Cone { " << endl;
	os << "                      bottomRadius  " << rr << endl;
	os << "                      height " << h << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;
}
