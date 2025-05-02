#include "mcEBTSourcePhsp.h"
#include "mcDefs.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "../geometry/text.h"

mcEBTSourcePhsp::mcEBTSourcePhsp(const char* name, int nThreads, double z)
	:mcSource(name, nThreads)
	, z_(z)
	, t_(mc_particle_t::MCP_PHOTON)
	, nparticles_(0)
	, idx_(nThreads, 0)
	, npars_(7)
{
}

void mcEBTSourcePhsp::sample(mcParticle& p, mcThread* thread)
{
	mcRng& rng = thread->rng();

	int i = ((idx_[thread->id()]++) % nparticles_) * npars_;
	p.t = t_;
	p.ke = phsp_[i++];
	p.weight = 1.0;
	p.thread_ = thread;

	double phi = rng.rnd() * TWOPI;
	double sina = sin(phi);
	double cosa = cos(phi);

	double x = phsp_[i++];
	double y = phsp_[i++];
	p.p.set(x * cosa - y * sina, x * sina + y * cosa, phsp_[i++] + z_);

	x = phsp_[i++];
	y = phsp_[i++];
	p.u.set(x * cosa - y * sina, x * sina + y * cosa, phsp_[i++]);

	etotal_[thread->id()] += p.ke * p.weight;
}

//void mcEBTSourcePhsp::loadData(istream& is)
//{
//	vector<double> data(npars_, 0);
//	string line;
//	std::getline(is, line, '\n');
//	phsp_.clear();
//
//	while (!is.fail())
//	{
//		std::getline(is, line, '\n');
//		if (is.fail() || line.length() < 30)
//			break;
//		if (GetFloatArray(line, data) != npars_)
//			throw std::exception("mcEBTSourcePhsp::loadData: wrong parameters number");
//		for (int i = 0; i < npars_; i++)
//			phsp_.push_back(data[i]);
//	}
//	nparticles_ = (int)phsp_.size() / npars_;
//}


void mcEBTSourcePhsp::readFromMemory(void* buffer)
{
	char* pbuffer = (char*)buffer;
	nparticles_ = *(int*)(pbuffer); pbuffer += sizeof(int);
	phsp_.resize(nparticles_ * npars_, 0);
	memcpy(&phsp_[0], pbuffer, nparticles_ * npars_ * sizeof(double));
}

void mcEBTSourcePhsp::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	double r = 0.025, rr = 0.1, h = 0.25, L = 0.5;

	os << "# Source: " << name_ << endl;

	geomVector3D p(0, 0, z_ - L * 0.5 - h);

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

	p(2) = z_ - h * 0.5;

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
