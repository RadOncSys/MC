#include "mcScoreTrack.h"
#include "mcGeometry.h"

mcScoreTrack::mcScoreTrack(int nThreads, double R, double Z1, double Z2, double EMIN, bool doPhotons, bool doElectrons, bool doPositrons) :
	mcScore(nullptr, nThreads), R_(R), Z1_(Z1), Z2_(Z2), EMIN_(EMIN),
	doPhotons_(doPhotons), doElectrons_(doElectrons), doPositrons_(doPositrons)
{
	photons_.resize(nThreads);
	electrons_.resize(nThreads);
	positrons_.resize(nThreads);
}

mcScoreTrack::~mcScoreTrack(void)
{
}

void mcScoreTrack::score(int iThread
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1
	, double ke)
{
	if ((pt == MCP_PHOTON && !doPhotons_) || (pt == MCP_NEGATRON && !doElectrons_) || (pt == MCP_POSITRON && !doPositrons_) ||
		ke < EMIN_)
		return;

	geomVector3D pp0 = p0;
	geomVector3D pp1 = p1;

	double r0 = pp0.lengthXY();
	double r1 = pp1.lengthXY();
	if (r0 >= R_ && r1 >= R_)
		return;

	if (r0 > R_)
	{
		geomVector3D v = pp0 - pp1;
		v.normalize();
		double d = mcGeometry::getDistanceToInfiniteCylinderInside(pp1, v, R_);
		pp0 = pp1 + (v * d);
	}
	else if (r1 > R_)
	{
		geomVector3D v = pp1 - pp0;
		v.normalize();
		double d = mcGeometry::getDistanceToInfiniteCylinderInside(pp0, v, R_);
		pp1 = pp0 + (v * d);
	}

	double z0 = pp0.z(), z1 = pp1.z();
	if ((z0 <= Z1_ && z1 <= Z1_) || (z0 >= Z2_ && z1 >= Z2_))
		return;
	if ((pp0 - pp1).sqLength() < 0.0001)
		return;

	if (z0 < Z1_)
	{
		geomVector3D v = pp0 - pp1;
		v.normalize();
		double d = (Z1_ - z1) / v.z();
		pp0 = pp1 + (v * d);
	}
	else if (z1 < Z1_)
	{
		geomVector3D v = pp1 - pp0;
		v.normalize();
		double d = (Z1_ - z0) / v.z();
		pp1 = pp0 + (v * d);
	}

	if (z0 > Z2_)
	{
		geomVector3D v = pp0 - pp1;
		v.normalize();
		double d = (Z2_ - z1) / v.z();
		pp0 = pp1 + (v * d);
	}
	else if (z1 > Z2_)
	{
		geomVector3D v = pp1 - pp0;
		v.normalize();
		double d = (Z2_ - z0) / v.z();
		pp1 = pp0 + (v * d);
	}

	if (pt == MCP_PHOTON) {
		photons_[iThread].push_back(pp0);
		photons_[iThread].push_back(pp1);
	}
	else if (pt == MCP_NEGATRON) {
		electrons_[iThread].push_back(pp0);
		electrons_[iThread].push_back(pp1);
	}
	else if (pt == MCP_POSITRON) {
		positrons_[iThread].push_back(pp0);
		positrons_[iThread].push_back(pp1);
	}
}

void mcScoreTrack::dumpVRML(ostream& os) const
{
	int i, j, ia;

	os << "# Particles tracks: " << this->getName() << endl;
	os << "Group {" << endl;
	os << "  children [" << endl;

	for (ia = 0; ia < 4; ia++)
	{
		const vector<vector<geomVector3D>>& t = ia == 0 ? photons_ : ia == 1 ? electrons_ : positrons_;
		if (t.empty())
			continue;

		os << "    Transform {" << endl;
		os << "      children Shape {" << endl;
		os << "        appearance Appearance {" << endl;
		os << "          material Material {" << endl;

		if (ia == 0)
			os << "            emissiveColor 0 0 1" << endl;
		else if (ia == 1)
			os << "            emissiveColor 0 1 0" << endl;
		else if (ia == 2)
			os << "            emissiveColor 0 0.5 1" << endl;
		else
			os << "            emissiveColor 1 0 0" << endl;

		//if(ia == 0)
		//	os << "            diffuseColor 0 1 0" << endl;
		//else if(ia == 1)
		//	os << "            diffuseColor 0 0 1" << endl;
		//else if(ia == 2)
		//	os << "            diffuseColor 0 0.5 1" << endl;
		//else
		//	os << "            diffuseColor 1 0 0" << endl;

		os << "          }" << endl;
		os << "        }" << endl;
		os << "        geometry IndexedLineSet {" << endl;
		os << "            coord Coordinate {" << endl;
		os << "                point [" << endl;

		for (j = 0; j < (int)t.size(); j++)
		{
			for (i = 0; i < (int)t[j].size(); i++)
			{
				os << "                    " << t[j][i].x() << ' ' << t[j][i].y() << ' ' << t[j][i].z() << ", " << endl;
				i++;
				os << "                    " << t[j][i].x() << ' ' << t[j][i].y() << ' ' << t[j][i].z();
				if (i < (int)t[j].size() - 1) os << ", ";
				os << endl;
			}
		}

		os << "                ]" << endl;
		os << "            }" << endl;
		os << "            coordIndex [" << endl;

		int jcount = 0;
		for (j = 0; j < (int)t.size(); j++)
		{
			for (i = 0; i < (int)t[j].size(); i += 2, jcount += 2)
			{
				os << "                " << jcount << " " << (jcount + 1) << " -1" << endl;
			}
		}

		os << "            ]" << endl;
		os << "        }" << endl;
		os << "      }" << endl;
		os << "    }" << endl;
	}

	os << "  ]" << endl;
	os << "}" << endl;
}
