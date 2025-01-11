#include "mcSourceProtonRF.h"
#include "mcDefs.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"
#include "../geometry/text.h"

mcSourceProtonRF::mcSourceProtonRF(const char* name, int nThreads, double z)
	:mcSource(name, nThreads)
	, z_(z)
	, nparticles_(0)
	, idx_(nThreads, 0)
{
	rng_.resize(nThreads);
	for (int i = 0; i < nThreads; i++)
		rng_[i].init(33, i + 1);
}

void mcSourceProtonRF::sample(mcParticle& p, mcThread* thread)
{
	int ithread = thread->id();
	int i = (idx_[ithread]++) % nparticles_;
	p = particles_[i];
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[ithread] += p.ke * p.weight;

	// Предполагаем радиальную симметрию и при каждом самплинге 
	// случайно поворачиваем частицу вокруг оси Z.

	double a = 2 * PI * rng_[ithread].rnd();
	double sa = sin(a), ca = cos(a);
	double x = p.p.p_[0], y = p.p.p_[1];
	p.p.p_[0] = x * ca + y * sa;
	p.p.p_[1] = -x * sa + y * ca;

	x = p.u.p_[0]; y = p.u.p_[1];
	p.u.p_[0] = x * ca + y * sa;
	p.u.p_[1] = -x * sa + y * ca;
}

void mcSourceProtonRF::loadData(istream& is)
{
	mcParticle particle;
	particle.t = mc_particle_t::MCP_PROTON;
	particle.q = 1;
	particle.thread_ = nullptr;
	particle.trackScore_ = nullptr;

	vector<double> data(5, 0);
	string line;

	// Первую строку пропускаем
	std::getline(is, line, '\n');

	while (!is.fail())
	{
		std::getline(is, line, '\n');
		if (is.fail() || line.length() < 30)
			break;
		if (GetFloatArray(line, data) != 5)
			throw std::exception("mcSourceProtonRF::loadData: 8 parameters expected");

		particle.ke = data[4];
		particle.p.set(data[0] * 0.1, data[2] * 0.1, z_);
		particle.u.set(data[1], data[3], sqrt( 1 - data[1] * data[1] - data[3] * data[3]));
		particle.weight = 1.0;

		particles_.push_back(particle);
	}
	nparticles_ = (int)particles_.size();
}

void mcSourceProtonRF::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	double r = 0.1, rr = 0.4, h = 1.0, L = 2.0;

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
