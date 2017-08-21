#include "mcSourceXrayBeam.h"
#include "mcDefs.h"
#include "mcSamplers.h"
#include "mcHistogramSampler.h"
#include "mcThread.h"
#include "../geometry/vec3d.h"

// Spectra
double ebins_06X[] = { 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5 };
double spec_06X[] = { 0.0186198, 0.0399641, 0.0405027, 0.0385809, 0.0362754, 0.050912, 0.058875, 0.0503812, 0.0634955, 0.0545269, 0.0264792 };

double ebins_10X[] = { 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9 };
double spec_10X[] = { 3.89224E-05, 0.013095, 0.0207605, 0.0213944, 0.0203683, 0.0285345, 0.0328886, 0.028523, 0.0372739, 0.0379239, 0.0287197, 0.0211778, 0.0147925, 0.00925699, 0.00437173 };

double ebins_18X[] = { 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16 };
double spec_18X[] = { 1.41302E-41, 4.97089E-07, 0.00109451, 0.00599288, 0.00979936, 0.0168469, 0.0216932, 0.0188831, 0.0242913, 0.0242057, 0.0186543, 0.0147144, 0.0117833, 0.00951739, 0.00771155, 0.00935528, 0.00793689, 0.00460063, 0.00203667};

mcSourceXrayBeam::mcSourceXrayBeam(const char* name, int nThreads, double z, double enom, double sad, double x1, double x2, double y1, double y2)
	: mcSource(name, nThreads)
	, esampler_(nullptr)
	, type_(MCP_PHOTON)
	, z_(z)
	, nom_energy_(enom)
	, sad_(sad)
	, fsx1_(x1)
	, fsx2_(x2)
	, fsy1_(y1)
	, fsy2_(y2)
{
	esampler_ = new mcHistogramSampler();

	if (nom_energy_ == 6.0)
		esampler_->setFromDistribution(100, sizeof(spec_06X) / sizeof(double), ebins_06X, spec_06X);
	else if (nom_energy_ == 10.0)
		esampler_->setFromDistribution(100, sizeof(spec_10X) / sizeof(double), ebins_10X, spec_10X);
	else if (nom_energy_ == 18.0)
		esampler_->setFromDistribution(100, sizeof(spec_18X) / sizeof(double), ebins_18X, spec_18X);
	else 
		throw std::exception("mcSourceXrayBeam: only 6, 10 and 18 MV theraputic photon beams supported");
}

mcSourceXrayBeam::~mcSourceXrayBeam(void)
{
	if (esampler_) delete esampler_;
}

void mcSourceXrayBeam::sample(mcParticle& p, mcThread* thread)
{
	// Согласование сигма по радиусу и отдельной координате.
	const double s = sqrt(0.5);
	mcRng& rng = thread->rng();

	p.t = type_;
	p.q = 0;
	p.ke = esampler_->sample(rng.rnd());

	double x = fsx1_ + (fsx2_ - fsx1_) * rng.rnd();
	double y = fsy1_ + (fsy2_ - fsy1_) * rng.rnd();

	p.p.set(0, 0, z_);
	p.u.set(x, y, sad_);
	p.u.normalize();

	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke;
}

void mcSourceXrayBeam::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	int i = 0;
	geomVector3D p[5];

	p[i++] = geomVector3D(0, 0, z_);
	p[i++] = geomVector3D(fsx1_, fsy1_, z_ + sad_);
	p[i++] = geomVector3D(fsx1_, fsy2_, z_ + sad_);
	p[i++] = geomVector3D(fsx2_, fsy2_, z_ + sad_);
	p[i++] = geomVector3D(fsx2_, fsy1_, z_ + sad_);

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

	for (i = 0; i < 5; i++) {
		os << "                    " << p[i].x() << ' ' << p[i].y() << ' ' << p[i].z();
		if (i < 4) os << ", ";
		os << endl;
	}

	os << "                ]" << endl;
	os << "            }" << endl;
	os << "            coordIndex [" << endl;

	os << "                0, 1, 2, -1," << endl;
	os << "                0, 2, 3, -1," << endl;
	os << "                0, 3, 4, -1," << endl;
	os << "                0, 4, 1, -1," << endl;

	os << "            ]" << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "    }" << endl;
}
