#include "mcClinicalElectronBeam.h"
#include "mcThread.h"
#include "mcDefs.h"
#include "mcGeometry.h"
#include "mcSamplers.h"
#include "../geometry/vec3d.h"

mcClinicalElectronBeam::mcClinicalElectronBeam(const char* name, int nThreads, double ke, double z0,
	double sigmaE, double sigmaA, double fsx, double fsy)
	: mcSource(name, nThreads)
	, type_(MCP_NEGATRON), ke_(ke), q_(-1)
	, sigmaE_(sigmaE), sigmaA_(sigmaA), fsx_(fsx), fsy_(fsy)
{
	p_.set(0, 0, z0);
	v_.set(0, 0, 1);
}

void mcClinicalElectronBeam::sample(mcParticle& p, mcThread* thread)
{
	const double scd = 0.95;	// масштаб перехода от изоцентра к рамке электронного аппликатора
	mcRng& rng = thread->rng();
	double x = (rng.rnd() - 0.5) * fsx_ * scd;
	double y = (rng.rnd() - 0.5) * fsy_ * scd;
	p.p.set(x, y, p_.z() - 5.0);

	geomVector3D v(x, y, scd * 100);
	v.normalize();

	double r = tan(sigmaA_ * mcSamplers::SampleGauss(rng.rnd()));
	double phi = rng.rnd() * TWOPI;
	double dx = r * cos(phi);
	double dy = r * sin(phi);
	v += geomVector3D(dx, dy, 0);
	v.normalize();
	p.u = v;

	double f = mcSamplers::SampleGauss(rng.rnd());
	if (fabs(f) > 2) f = 0;	// Таких частиц мало. Просто подменяем их монохроматичными, чтобы не усложнять код
	p.ke = ke_ + 0.5 * sigmaE_ * f;

	p.t = type_;
	p.q = q_;
	p.plast = p.p;
	p.weight = 1;
	p.thread_ = thread;
	p.trackScore_ = trackScore_;
	etotal_[thread->id()] += p.ke;
}

void mcClinicalElectronBeam::dumpVRML(ostream& os) const
{
	mcSource::dumpVRML(os);

	double r0 = 0.4, r = 0.2, rr = 0.4, hh = 1.0, L = 5.0;
	double sad = 100, cad = 5;

	os << "# Source: " << name_ << endl;

	// Шарик с упирающейся стрелкой

	os << "Transform {" << endl;
	os << "  translation " << p_.x() << ' ' << p_.y() << ' ' << (p_.z() - sad) << endl;
	os << "  children [" << endl;
	os << "    Shape{" << endl;
	os << "      appearance Appearance {" << endl;
	os << "        material Material {" << endl;
	os << "          diffuseColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "          transparency " << transparancy_ << endl;
	os << "        }" << endl;
	os << "      }" << endl;
	os << "      geometry Sphere { radius " << r0 << " }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;

	geomVector3D p = p_ + (v_ * (L*0.5));

	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << (p_.z() - sad) << endl;
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

	p = p_ + (v_ * (L + hh * 0.5));

	os << "Transform {" << endl;
	os << "  translation " << p.x() << ' ' << p.y() << ' ' << (p_.z() - sad) << endl;
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
	os << "                      height " << hh << endl;
	os << "      }" << endl;
	os << "    }" << endl;
	os << "  ]" << endl;
	os << "}" << endl;

	// Рамка электронного аппликатора

	os << "# Rectangle ring: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	double d = 5, h = 1;
	double x1 = -fsx_ / 2, x2 = -x1;
	double y1 = -fsy_ / 2, y2 = -y1;

	geomVector3D pp[16];
	pp[0] = geomVector3D(x1 - d, y1 - d, -h-cad);
	pp[1] = geomVector3D(x1 - d, y2 + d, -h-cad);
	pp[2] = geomVector3D(x2 + d, y2 + d, -h-cad);
	pp[3] = geomVector3D(x2 + d, y1 - d, -h-cad);
	pp[4] = geomVector3D(x1, y1, -h-cad);
	pp[5] = geomVector3D(x1, y2, -h-cad);
	pp[6] = geomVector3D(x2, y2, -h-cad);
	pp[7] = geomVector3D(x2, y1, -h-cad);
	pp[8] = geomVector3D(x1 - d, y1 - d, - cad);
	pp[9] = geomVector3D(x1 - d, y2 + d, - cad);
	pp[10] = geomVector3D(x2 + d, y2 + d, - cad);
	pp[11] = geomVector3D(x2 + d, y1 - d, - cad);
	pp[12] = geomVector3D(x1, y1, -cad);
	pp[13] = geomVector3D(x1, y2, -cad);
	pp[14] = geomVector3D(x2, y2, -cad);
	pp[15] = geomVector3D(x2, y1, -cad);

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
		os << "                    " << pp[i].x() << ' ' << pp[i].y() << ' ' << pp[i].z();
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

	os << "  ]" << endl;
	os << "}" << endl;
}
