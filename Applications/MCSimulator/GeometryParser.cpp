#include "GeometryParser.h"
#include "XmlParseReaderBase.h"

#include "../../mc/mc/mcTransportCylinder.h"
#include "../../mc/mc/mcTransportCone.h"
#include "../../mc/mc/mcTransportPrism.h"
#include "../../mc/mc/mcTransportRing.h"
#include "../../mc/mc/mcTransportConicalRing.h"
#include "../../mc/mc/mcTransportConicalHole.h"
#include "../../mc/mc/mcTransportRectangleRing.h"
#include "../../mc/mc/mcTransportJawPairRounded.h"
#include "../../mc/mc/mcTransportRectanglePolygonSideHole.h"
#include "../../mc/mc/mcTransportMLC.h"
#include "../../mc/mc/mcTransportCylinderStack.h"
#include "../../mc/mc/mcTransportSlab.h"
#include "../../mc/mc/mcTransportPlaneFilter.h"
#include "../../mc/mc/mcTransportRectangleTrap.h"
#include "../../mc/mc/mcTransportWedge.h"
#include "../../mc/mc/mcTransportAxialSymmetricSplitter.h"
#include "../../mc/mc/mcTransportSimpleSplitter.h"
#include "../../mc/mc/mcPTLasVegas.h"
#include "../../mc/mc/mcETransportSphere.h"
#include "../../mc/mc/mcETransportTrap.h"
#include "../../mc/mc/mcETransportConvexPolygonCircle.h"
#include "../../mc/mc/mcTransportEmbeddedGroup.h"
#include "../../mc/mc/mcTransportLinearChain.h"

#include "../../mc/mc/mcScorePHSP.h"
#include "../../mc/mc/mcScoreBeamFluence.h"
#include "../../mc/mc/mcScoreBeamFluence2.h"
#include "../../mc/mc/mcScoreMatrixRZ.h"
#include "../../mc/mc/mcScoreConicalRZ.h"
#include "../../mc/mc/mcScoreBeamFluenceXY.h"
#include "../../mc/mc/mcScorePhaseSpaceConcentrator.h"
#include "../../mc/mc/mcScoreMatrixXY.h"
#include "../../mc/mc/mcScoreMatrix2D.h"
#include "../../mc/mc/mcScoreEnergySpectrum.h"
#include "../../mc/mc/mcScoreEnergyFluence.h"
#include "../../mc/mc/mcScoreSpectraFluence.h"
#include "../../mc/mc/mcScoreParticleContainer.h"
#include "../../mc/mc/mcScoreSphereFluence.h"
#include "../../mc/mc/mcScoreSphereMatrix.h"

#include "../../mc/mc/mcSourceSimpleMono.h"
#include "../../mc/mc/mcSourceAccelerator.h"
#include "../../mc/mc/mcSourceXraySigmaR.h"
#include "../../mc/mc/mcSourceCylindricalC60.h"
#include "../../mc/mc/mcSourceModelRadialPhsp.h"
#include "../../mc/mc/mcSourceLEBA.h"
#include "../../mc/mc/mcSourceSphereC60.h"
#include "../../mc/mc/mcSourceAcceleratedBeam.h"
#include "../../mc/mc/mcClinicalElectronBeam.h"

#include <io.h>
#include <fcntl.h>
#include <algorithm>
#include <fstream>

using namespace std;

enum mc_particle_t GeometryParser::convert_S2T_ptype(const wchar_t* st)
{
	if (_wcsicmp(st, L"photon") == 0 || _wcsicmp(st, L"c60") == 0 || _wcsicmp(st, L"c60_sphere") == 0 || _wcsicmp(st, L"phsp_photon") == 0)
		return MCP_PHOTON;
	else if (_wcsicmp(st, L"electron") == 0)
		return MCP_NEGATRON;
	else if (_wcsicmp(st, L"positron") == 0)
		return MCP_POSITRON;
	else
		throw exception("Unknown radiation type");
}

enum spectrum_distr_t GeometryParser::convert_spec_type(const wchar_t* st)
{
	if (_wcsicmp(st, L"gauss") == 0)
		return SPECTRUM_GAUSS;
	else if (_wcsicmp(st, L"triangle") == 0)
		return SPECTRUM_TRIANGLE;
	else if (_wcsicmp(st, L"prism") == 0)
		return SPECTRUM_PRISM;
	else
		throw exception("Unknown spectrum type");
}

enum profile_distr_t GeometryParser::convert_profile_distr_type(const wchar_t* st)
{
	if (_wcsicmp(st, L"gauss") == 0)
		return PROFILE_GAUSS;
	else if (_wcsicmp(st, L"exponent") == 0)
		return PROFILE_EXPONENT;
	else if (_wcsicmp(st, L"prism") == 0)
		return PROFILE_PRISM;
	else
		throw exception("Unknown beam profile type");
}

mcTransport* GeometryParser::ParseTransport(const XPRNode& geometry, const mcMedia* media, int nThreads)
{
	// Используется для обработки ситуации, когда модуль является оболочкой группы а не реальным модулем. 
	// В этом случае логика организована так, что функцию группового модуля выполняет один, 
	// который иницилизируется в своем цикле парсинга.
	bool skipInit = false;

	mcTransport* t = nullptr;

	// Общие для всех модулей параметры, которые должны быть получены в результате парсинга
	wstring geomType;
	wstring geomName;
	wstring geomMedium;
	double densityRatio = 1.0;
	double electron_cut = 0, photon_cut = 0;
	double cr = 1.0, cg = 1.0, cb = 1.0, ct = 1.0;
	double x0, y0, z0;
	double vx, vy, vz;
	double xx, xy, xz;

	// Специфичные для модуля параметры.
	// Для простоты перечисляем все возможные варианты.
	double radius = 0;
	double width = 0;
	double height = 0;
	double length = 0;
	double ax = 0, ay = 0, az = 0;
	double r0 = 0, r1 = 0;
	double dr = 0, ds = 0;
	double sx0 = 0, sy0 = 0;
	double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
	double d = 0;
	double h = 0;
	double focus = 0;
	int nc = 0;
	int nsplit = 0;
	mc_particle_t ptype = mc_particle_t::MCP_PHOTON;
	
	std::vector<double> poly_z;
	std::vector<double> poly_x;
	std::vector<double> poly_y;

	for (auto node : geometry.Nodes)
	{
		if (_wcsicmp(node.Name.c_str(), L"type") == 0)
			geomType = node.Text;
		else if (_wcsicmp(node.Name.c_str(), L"name") == 0)
			geomName = node.Text;
		else if (_wcsicmp(node.Name.c_str(), L"medium") == 0)
			geomMedium = node.Text;
		else if (_wcsicmp(node.Name.c_str(), L"density") == 0)
			densityRatio = _wtof(node.Text.c_str());

		// Transport Cutoff
		else if (_wcsicmp(node.Name.c_str(), L"transCutoff") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"electron") == 0)
					electron_cut = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"photon") == 0)
					photon_cut = _wtof(n1.Text.c_str());
			}
		}

		// Color
		else if (_wcsicmp(node.Name.c_str(), L"Color") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"r") == 0)
					cr = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"g") == 0)
					cg = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"b") == 0)
					cb = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"t") == 0)
					ct = _wtof(n1.Text.c_str());
			}
		}

		// Location
		else if (_wcsicmp(node.Name.c_str(), L"position") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"x") == 0)
					x0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y") == 0)
					y0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z") == 0)
					z0 = _wtof(n1.Text.c_str());
			}
		}
		else if (_wcsicmp(node.Name.c_str(), L"normal") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"x") == 0)
					vx = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y") == 0)
					vy = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z") == 0)
					vz = _wtof(n1.Text.c_str());
			}
		}
		else if (_wcsicmp(node.Name.c_str(), L"xaxis") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"x") == 0)
					xx = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y") == 0)
					xy = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z") == 0)
					xz = _wtof(n1.Text.c_str());
			}
		}

		// Специфичные для модуля параметры
		else if (_wcsicmp(node.Name.c_str(), L"size") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"radius") == 0)
					radius = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"width") == 0)
					width = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"height") == 0)
					height = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"length") == 0)
					length = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ax") == 0)
					ax = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ay") == 0)
					ay = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"az") == 0)
					az = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"r0") == 0)
					r0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"r1") == 0)
					r1 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"focus") == 0)
					focus = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"x0") == 0)
					sx0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y0") == 0)
					sy0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"x1") == 0)
					x1 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"x2") == 0)
					x2 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y1") == 0)
					y1 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y2") == 0)
					y2 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"d") == 0)
					d = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"h") == 0)
					h = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"nc") == 0)
					nc = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"dr") == 0)
					dr = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ds") == 0)
					ds = _wtof(n1.Text.c_str());
			}
		}

		else if (_wcsicmp(node.Name.c_str(), L"polygon") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"point") == 0)
				{
					double x = 0, z = 0;
					for (auto n2 : n1.Nodes)
					{
						if (_wcsicmp(n2.Name.c_str(), L"x") == 0)
							x = _wtof(n2.Text.c_str());
						else if (_wcsicmp(n2.Name.c_str(), L"z") == 0)
							z = _wtof(n2.Text.c_str());
					}
					poly_z.push_back(z);
					poly_x.push_back(x);
					poly_y.push_back(x);
				}
			}
		}

		else if (_wcsicmp(node.Name.c_str(), L"nsplit") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"particle") == 0)
					ptype = convert_S2T_ptype(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"value") == 0)
					nsplit = _wtoi(n1.Text.c_str());
			}
		}
	}

	//
	// Создание объектов
	//
	geomVector3D origin(x0, y0, z0);
	geomVector3D normal(vx, vy, vz);
	geomVector3D xaxis(xx, xy, xz);
	normal.normalize();
	xaxis.normalize();

	if (_wcsicmp(geomType.c_str(), L"cylinder") == 0)
	{
		t = new mcTransportCylinder(origin, normal, xaxis, radius, height);
	}
	else if (_wcsicmp(geomType.c_str(), L"cone") == 0)
	{
		t = new mcTransportCone(origin, normal, xaxis, radius, height);
	}
	else if (_wcsicmp(geomType.c_str(), L"prism") == 0)
	{
		t = new mcTransportPrism(origin, normal, xaxis, ax, ay, az);
	}
	else if (_wcsicmp(geomType.c_str(), L"wedge") == 0)
	{
		t = new mcTransportWedge(origin, normal, xaxis, ax, ay, az);
	}
	else if (_wcsicmp(geomType.c_str(), L"ring") == 0)
	{
		t = new mcTransportRing(origin, normal, xaxis, r0, r1, height);
	}
	else if (_wcsicmp(geomType.c_str(), L"conicalring") == 0)
	{
		t = new mcTransportConicalRing(origin, normal, xaxis, r0, r1, height, focus);
	}
	else if (_wcsicmp(geomType.c_str(), L"conicalhole") == 0)
	{
		t = new mcTransportConicalHole(origin, normal, xaxis, r0, r1, height, focus);
	}
	else if (_wcsicmp(geomType.c_str(), L"rectanglering") == 0)
	{
		t = new mcTransportRectangleRing(origin, normal, xaxis, x1, x2, y1, y2, d, h);
	}
	else if (_wcsicmp(geomType.c_str(), L"jaw_pair") == 0)
	{
		t = new mcTransportJawPairRounded(origin, normal, xaxis, radius, height, width);
		((mcTransportJawPairRounded*)t)->setFS(x1, x2);
	}
	else if (_wcsicmp(geomType.c_str(), L"rectanglepolygonsidehole") == 0)
	{
		std::for_each(poly_x.begin(), poly_x.end(), [sx0](double& x)
		{
			x += sx0;
		});
		std::for_each(poly_y.begin(), poly_y.end(), [sy0](double& y) { y += sy0; });
		t = new mcTransportRectanglePolygonSideHole(origin, normal, xaxis,
			d + sx0, d + sy0, poly_z, poly_x, poly_y);
	}
	else if (_wcsicmp(geomType.c_str(), L"mlc") == 0)
	{
		t = new mcTransportMLC(origin, normal, xaxis, focus, radius, height, width, length, x1, x2, y1, y2);
	}
	else if (_wcsicmp(geomType.c_str(), L"gridcylinder") == 0)
	{
		t = new mcTransportCylinderStack(origin, normal, xaxis, nc, dr, ds, height);
	}
	else if (_wcsicmp(geomType.c_str(), L"planefilter") == 0)
	{
		t = new mcTransportPlaneFilter(origin, normal, xaxis);
	}
	else if (_wcsicmp(geomType.c_str(), L"etrap") == 0)
	{
		t = new mcETransportTrap(origin, normal, xaxis);
	}
	else if (_wcsicmp(geomType.c_str(), L"rectangletrap") == 0)
	{
		t = new mcTransportRectangleTrap(origin, normal, xaxis, width, length);
	}
	else if (_wcsicmp(geomType.c_str(), L"axial_splitter") == 0)
	{
		t = new mcTransportAxialSymmetricSplitter(origin, normal, xaxis, ptype, nsplit);
	}
	else if (_wcsicmp(geomType.c_str(), L"simple_splitter") == 0)
	{
		t = new mcTransportSimpleSplitter(origin, normal, xaxis, ptype, nsplit);
	}
	else if (_wcsicmp(geomType.c_str(), L"slab") == 0)
	{
		t = new mcTransportSlab(origin, normal, xaxis, height);
	}
	else if (_wcsicmp(geomType.c_str(), L"group") == 0)
	{
		t = new mcTransport(origin, normal, xaxis);

		// Ищем и парсим вложенные модули
		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* region = GeometryParser::ParseTransport(node, media, nThreads);
				region->setDefaultScore(nThreads);
				t->addRegion(region);
			}
		}
	}
	else if (_wcsicmp(geomType.c_str(), L"embedding") == 0)
	{
		// Временный транспорт нужен для системы координат группы.
		mcTransport ttmp(origin, normal, xaxis);
		mcTransport* eprev = nullptr;

		// Ищем и парсим вложенные модули
		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* eobj = GeometryParser::ParseTransport(node, media, nThreads);
				if (eobj == nullptr)
					throw exception("Embedding group should contain only embedding objects");

				// Корректируем систему координат из группы в сцену.
				eobj->MoveToCoordinateSystem(ttmp.MT2W());

				eobj->setDefaultScore(nThreads);
				if (eprev != nullptr)
				{
					eprev->setExternalTransport(eobj);
					eobj->setInternalTransport(eprev);
				}
				eprev = eobj;
			}
		}
		t = eprev;
		skipInit = true;
	}
	else if (_wcsicmp(geomType.c_str(), L"embedded_group") == 0)
	{
		t = new mcTransportEmbeddedGroup(origin, normal, xaxis);

		// Ищем и парсим вложенные модули
		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* eobj = GeometryParser::ParseTransport(node, media, nThreads);
				if (eobj == nullptr)
					throw exception("Embedding group internal module parse error");

				// Корректируем систему координат из группы в сцену.
				eobj->MoveToCoordinateSystem(t->MT2W());
				((mcTransportEmbeddedGroup*)t)->addTransport(eobj);
			}
		}
	}
	else if (_wcsicmp(geomType.c_str(), L"linear_chain") == 0)
	{
		t = new mcTransportLinearChain(origin, normal, xaxis);

		// Ищем и парсим вложенные модули
		for (auto node : geometry.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* eobj = GeometryParser::ParseTransport(node, media, nThreads);
				if (eobj == nullptr)
					throw exception("Embedding group should contain only embedding objects");

				// Корректируем систему координат из группы в сцену.
				eobj->MoveToCoordinateSystem(t->MT2W());

				((mcTransportLinearChain*)t)->addTransport(eobj);
			}
		}
		((mcTransportLinearChain*)t)->completeInit();
	}
	else if (_wcsicmp(geomType.c_str(), L"ptlasvegas") == 0)
	{
		t = new mcPTLasVegas(origin, normal, xaxis);
	}
	else if (_wcsicmp(geomType.c_str(), L"esphere") == 0)
	{
		t = new mcETransportSphere(origin, normal, xaxis, radius);
	}
	else if (_wcsicmp(geomType.c_str(), L"e_convex_polygon_circle") == 0)
	{
		t = new mcETransportConvexPolygonCircle(origin, normal, xaxis, poly_z, poly_x);
	}

	if (t == nullptr)
		throw exception(XmlParseReaderBase::copyWStringToStlString(
			(wstring(L"Unsupported geometry type: ") + geomType).c_str()).c_str());

	if (!skipInit)
	{
		t->transCutoff_elec = electron_cut;
		t->transCutoff_phot = photon_cut;

		t->setMediaRef(media);
		t->setName(XmlParseReaderBase::copyWStringToStlString(geomName.c_str()).c_str());
		if (!geomMedium.empty())
			t->setMediumMono(XmlParseReaderBase::copyWStringToStlString(geomMedium.c_str()).c_str());
		t->setMediumMono(densityRatio);
		t->setColor(cr, cg, cb, ct);
	}

	return t;
}

mcScore* GeometryParser::ParseScore(const XPRNode& item, int nThreads)
{
	mcScore* score = nullptr;

	// Общие для всех модулей параметры, которые должны быть получены в результате парсинга
	wstring scoreType;
	string scoreModule;

	// Специфичные для модуля параметры.
	// Для простоты перечисляем все возможные варианты.
	string outfile;
	int nr = 0, nx = 0, ny = 0, nz = 0;
	int ne = 0, nebins = 0, nt = 0, na = 0;
	int np = 0, nm = 0;
	double rmax = 0, zmin = 0, zmax = 0;
	double r0 = 0, r1 = 0;
	double ziso = 0, sad = 0;
	double emax = 0;
	double ecut = 0;
	int nr_s = 0;
	double rmax_s = 0, tmax = 0;
	double H = 0;
	double d_iz = 0;
	double psx = 0, psy = 0, psz = 0;
	double z0 = 0, focus = 0;
	double x1 = 0, y1 = 0, z1 = 0;
	double x2 = 0, y2 = 0, z2 = 0;
	bool isxray = false;
	mc_particle_t pt = mc_particle_t::MCP_PHOTON;
	wstring outdcm;
	wstring calibrdcm;

	// Структура для накопления данных для последующего создания подисточников
	struct ss_data
	{
		mc_particle_t pt;
		int na;
		double zmin, zmax, focus, amax;
		string ssname;
	};
	vector<ss_data> subsources;

	// Контрольные флаги
	bool isSizeSetFound = false;
	bool isSpecParsDefined = false;

	for (auto node : item.Nodes)
	{
		if (_wcsicmp(node.Name.c_str(), L"type") == 0)
			scoreType = node.Text;
		else if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			scoreModule = XmlParseReaderBase::copyWStringToStlString(node.Text.c_str());
		else if (_wcsicmp(node.Name.c_str(), L"outfile") == 0)
			outfile = XmlParseReaderBase::copyWStringToStlString(node.Text.c_str());
		else if (_wcsicmp(node.Name.c_str(), L"size") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"nr") == 0)
					nr = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"nx") == 0)
					nx = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ny") == 0)
					ny = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"nz") == 0)
					nz = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ne") == 0)
					ne = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"nebins") == 0)
					nebins = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"nt") == 0)
					nt = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"na") == 0)
					na = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"np") == 0)
					np = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"nm") == 0)
					nm = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"r0") == 0)
					r0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"r1") == 0)
					r1 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"rmax") == 0)
					rmax = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"zmin") == 0)
					zmin = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"zmax") == 0)
					zmax = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ziso") == 0)
					ziso = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"sad") == 0)
					sad = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"pt") == 0)
					pt = GeometryParser::convert_S2T_ptype(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"emax") == 0)
					emax = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"nr_s") == 0)
					nr_s = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"rmax_s") == 0)
					rmax_s = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"tmax") == 0)
					tmax = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"H") == 0)
					H = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"d_iz") == 0)
					d_iz = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"psx") == 0)
					psx = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"psy") == 0)
					psy = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"psz") == 0)
					psz = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z0") == 0)
					z0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"x1") == 0)
					x1 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y1") == 0)
					y1 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z1") == 0)
					z1 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"x2") == 0)
					x2 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y2") == 0)
					y2 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z2") == 0)
					z2 = _wtof(n1.Text.c_str());
			}
			if (node.Nodes.size() > 0)
				isSizeSetFound = true;
		}
		else if (_wcsicmp(node.Name.c_str(), L"spectrum") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"nebins") == 0)
					nebins = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"emax") == 0)
					emax = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ecut") == 0)
					ecut = _wtoi(n1.Text.c_str());
			}
			if (node.Nodes.size() > 0)
				isSpecParsDefined = true;
		}
		else if (_wcsicmp(node.Name.c_str(), L"image") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"outdcm") == 0)
					outdcm = n1.Text;
				else if (_wcsicmp(n1.Name.c_str(), L"calibrdcm") == 0)
					calibrdcm = n1.Text;
			}
		}
		else if (_wcsicmp(node.Name.c_str(), L"subsource") == 0)
		{
			ss_data ssrc;
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"name") == 0)
					ssrc.ssname = XmlParseReaderBase::copyWStringToStlString(n1.Name.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"amax") == 0)
					ssrc.amax = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"na") == 0)
					ssrc.na = _wtoi(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"zmin") == 0)
					ssrc.zmin = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"zmax") == 0)
					ssrc.zmax = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"focus") == 0)
					ssrc.focus = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"particle") == 0)
					ssrc.pt = GeometryParser::convert_S2T_ptype(n1.Text.c_str());
			}
			subsources.push_back(ssrc);
		}
		else if (_wcsicmp(node.Name.c_str(), L"isxray") == 0)
			isxray = _wcsicmp(node.Text.c_str(), L"true") == 0 ? true : false;

		else if (_wcsicmp(node.Name.c_str(), L"pt") == 0)
			pt = GeometryParser::convert_S2T_ptype(node.Text.c_str());
		else if (_wcsicmp(node.Name.c_str(), L"particles") == 0)
			pt = (mc_particle_t) _wtoi(node.Text.c_str());
		else if (_wcsicmp(node.Name.c_str(), L"focus") == 0)
			focus = _wtof(node.Text.c_str());
	}

	if (!isSizeSetFound && _wcsicmp(scoreType.c_str(), L"fluence_sphere") != 0 &&
		_wcsicmp(scoreType.c_str(), L"fluence_plane") != 0)
		throw exception(("Please, indicate score size for module: " + scoreModule).c_str());

	if (_wcsicmp(scoreType.c_str(), L"phsp") == 0)
	{
		if (outfile.empty())
			throw exception("Please, specify dump file for scoring");
		score = new mcScorePHSP(scoreModule.c_str(), outfile.c_str());
	}

	// 3D дозовое распределение в RZ геометрии
	else if (_wcsicmp(scoreType.c_str(), L"rz") == 0)
	{
		score = new mcScoreMatrixRZ(scoreModule.c_str(), nThreads, nr, nz, rmax, zmin, zmax);
	}

	// 3D дозовое распределение в веерной RZ геометрии
	else if (_wcsicmp(scoreType.c_str(), L"rz_conical") == 0)
	{
		score = new mcScoreConicalRZ(scoreModule.c_str(), nThreads, nr, nz, rmax, zmin, zmax, ziso, sad);
	}

	// Распределение потока и спектров в плоскости через кольца
	else if (_wcsicmp(scoreType.c_str(), L"e_spectrum") == 0)
	{
		score = new mcScoreEnergySpectrum(scoreModule.c_str(), nThreads, pt, ne, emax, rmax);
	}

	// Распределение потока энергии в плоскости через кольца
	else if (_wcsicmp(scoreType.c_str(), L"e_fluence") == 0)
	{
		score = new mcScoreEnergyFluence(scoreModule.c_str(), nThreads, pt, nr, rmax);
	}

	// Распределение потока энергии в плоскости через кольца
	else if (_wcsicmp(scoreType.c_str(), L"spectra_fluence") == 0)
	{
		score = new mcScoreSpectraFluence(scoreModule.c_str(), nThreads, pt, ecut, nr, nebins, rmax, emax);
	}

	// Поток излучения в RZ геометрии
	else if (_wcsicmp(scoreType.c_str(), L"fluence") == 0)
	{
		score = new mcScoreBeamFluence(scoreModule.c_str(),	nThreads, nr, rmax);
		for (auto ss : subsources)
		{
			((mcScoreBeamFluence*) score)->addSubsource(new mcScoreBeamFluenceSubsource(
				ss.ssname.c_str(), nThreads, ss.pt, ss.zmin, ss.zmax, ss.focus, nr, rmax, ne, emax, ss.na, ss.amax));
		}
	}
	else if (_wcsicmp(scoreType.c_str(), L"fluence2") == 0)
	{
		score = new mcScoreBeamFluence2(scoreModule.c_str(), nThreads, nr, rmax, nr_s, rmax_s, H, nz, d_iz);
	}

	// Дозовое распределение в плоскости XY
	else if (_wcsicmp(scoreType.c_str(), L"xy") == 0)
	{
		mcScoreMatrixXY* scorexy = new mcScoreMatrixXY(scoreModule.c_str(), nThreads, nx, ny, psx, psy, psz, z0);

		// Файлы изображений DICOM
		if (!outdcm.empty())
			scorexy->SetImageFile(XmlParseReaderBase::copyWStringToStlString(outdcm.c_str()).c_str());
		if (!calibrdcm.empty())
			scorexy->SetCalibrationFile(XmlParseReaderBase::copyWStringToStlString(calibrdcm.c_str()).c_str());
		score = scorexy;
	}

	else if (_wcsicmp(scoreType.c_str(), L"2D") == 0)
	{
		score = new mcScoreMatrix2D(scoreModule.c_str(), nThreads, nx, nz, geomVector3D(x1, y1, z1), geomVector3D(x2, y2, z2));
	}

	// Поток излучения в координатах XY
	else if (_wcsicmp(scoreType.c_str(), L"fluence_xy") == 0)
	{
		if (!isSpecParsDefined)
			throw exception("Please, indicate spectrum parameters");
		score = new mcScoreBeamFluenceXY(scoreModule.c_str(), nThreads, nx, ny, psx, psy, nebins, emax);
	}

	// 4-x мерное (энергия / радиус / отклонение / азимут) распределение частиц в плоскости
	else if (_wcsicmp(scoreType.c_str(), L"phsp_concentrator") == 0)
	{
		if (outfile.empty())
			throw exception("Please, specify file name for scoring");
		score = new mcScorePhaseSpaceConcentrator(scoreModule.c_str(), nThreads,
			pt, isxray, focus, ne, nr, nt, na, emax, rmax, tmax, outfile.c_str());
	}

	else if (_wcsicmp(scoreType.c_str(), L"fluence_plane") == 0)
	{
		auto sc = new mcScoreParticleContainer(scoreModule.c_str(), nThreads);
		sc->setParticleTypeFilter(pt);
		score = sc;
	}

	// простейший накопитель частиц вылетающих из сферы
	// (на самом деле может использоваться как dump всех фотонов вылетающих 
	// из последнего объекта перед ловушкой и зарегистрированных в момент вылета)
	else if (_wcsicmp(scoreType.c_str(), L"fluence_sphere") == 0)
	{
		score = new mcScoreSphereFluence(scoreModule.c_str(), nThreads);
	}

	else if (_wcsicmp(scoreType.c_str(), L"matrix_sphere") == 0)
	{
		score = new mcScoreSphereMatrix(scoreModule.c_str(), nThreads, np, nm, r0, r1);
	}

	else
		throw exception(XmlParseReaderBase::copyWStringToStlString((L"Unsupported score type: " + scoreType).c_str()).c_str());

	score->setColor(0.0, 0.2, 0.0);
	score->setName(XmlParseReaderBase::copyWStringToStlString(scoreType.c_str()).c_str());

	return score;
}

mcSource* GeometryParser::ParseSource(const XPRNode& item, int nThreads)
{
	mcSource* source = nullptr;

	// Общие для всех модулей параметры, которые должны быть получены в результате парсинга
	string srcName;
	string sourceModule;
	wstring radTypeName;
	double energy = 0, ewidth = 0, awidth = 0;
	spectrum_distr_t sptype = spectrum_distr_t::SPECTRUM_GAUSS;
	profile_distr_t prf_type = profile_distr_t::PROFILE_EXPONENT;
	double x0 = 0, y0 = 0, z0 = 0;
	double vx = 0, vy = 0, vz = 0;
	bool isStartInside = false;

	// Специфичные для модуля параметры.
	// Для простоты перечисляем все возможные варианты.
	wstring fileName;
	double d = 0, h = 0, radius = 0, size = 0, angle = 0;
	double fsx = 0, fsy = 0;
	wstring srcType;
	wstring shapeType;

	// Контрольные флаги
	bool isRadiationDefined = false;
	bool isPositionDefined = false;
	bool isDirectionDefined = false;
	bool isCoreSetDefined = false;
	bool isShapeDefined = false;

	for (auto node : item.Nodes)
	{
		if (_wcsicmp(node.Name.c_str(), L"name") == 0)
			srcName = XmlParseReaderBase::copyWStringToStlString(node.Text.c_str());
		else if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			sourceModule = XmlParseReaderBase::copyWStringToStlString(node.Text.c_str());
		else if (_wcsicmp(node.Name.c_str(), L"IsStartInside") == 0)
			isStartInside = _wcsicmp(node.Text.c_str(), L"true") == 0;

		// Radiation
		else if (_wcsicmp(node.Name.c_str(), L"radiation") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"type") == 0)
					radTypeName = n1.Text;
				else if (_wcsicmp(n1.Name.c_str(), L"energy") == 0)
					energy = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"spectrum") == 0)
					sptype = GeometryParser::convert_spec_type(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"ewidth") == 0)
					ewidth = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"awidth") == 0)
					awidth = _wtof(n1.Text.c_str());
			}
			if (node.Nodes.size() > 0)
				isRadiationDefined = true;
		}

		// Location
		else if (_wcsicmp(node.Name.c_str(), L"position") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"x") == 0)
					x0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y") == 0)
					y0 = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z") == 0)
					z0 = _wtof(n1.Text.c_str());
			}
			if (node.Nodes.size() > 0)
				isPositionDefined = true;
		}

		// Direction
		else if (_wcsicmp(node.Name.c_str(), L"direction") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"x") == 0)
					vx = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"y") == 0)
					vy = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"z") == 0)
					vz = _wtof(n1.Text.c_str());
			}
			if (node.Nodes.size() > 0)
				isDirectionDefined = true;
		}

		else if (_wcsicmp(node.Name.c_str(), L"phsp") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"file") == 0)
					fileName = n1.Text.c_str();
			}
		}

		else if (_wcsicmp(node.Name.c_str(), L"core") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"d") == 0)
					d = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"h") == 0)
					h = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"radius") == 0)
					radius = _wtof(n1.Text.c_str());
			}
			if (node.Nodes.size() > 0)
				isCoreSetDefined = true;
		}

		else if (_wcsicmp(node.Name.c_str(), L"shape") == 0)
		{
			for (auto n1 : node.Nodes)
			{
				if (_wcsicmp(n1.Name.c_str(), L"direction") == 0)
					srcType = n1.Text;
				else if (_wcsicmp(n1.Name.c_str(), L"shapetype") == 0)
					shapeType = n1.Text;
				else if (_wcsicmp(n1.Name.c_str(), L"size") == 0)
					size = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"angle") == 0)
					angle = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"profile") == 0)
					prf_type = GeometryParser::convert_profile_distr_type(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"fsx") == 0)
					fsx = _wtof(n1.Text.c_str());
				else if (_wcsicmp(n1.Name.c_str(), L"fsy") == 0)
					fsy = _wtof(n1.Text.c_str());
			}
			if (node.Nodes.size() > 0)
				isShapeDefined = true;
		}
	}

	if (sourceModule.empty())
		throw exception("Please, specify attached module for source");

	if (!isRadiationDefined)
		throw exception("Please, indicate radiation type");

	if (!isPositionDefined)
		throw exception("Please, indicate radiation source position");

	if (!isDirectionDefined)
		throw exception("Please, indicate radiation source direction");

	if (_wcsicmp(radTypeName.c_str(), L"phsp_photon") == 0)
	{
		// Источник на основе модели фазового пространства (гистограм распределений частиц)
		if (!isRadiationDefined)
			throw exception("Please, indicate phase space data");

		if (fileName.empty())
			throw exception("Please, indicate phase space model file name");

		// Read model data file
		int ifile;
		if (_sopen_s(&ifile, XmlParseReaderBase::copyWStringToStlString(fileName.c_str()).c_str(),
			_O_RDONLY | _O_SEQUENTIAL | _O_BINARY, _SH_DENYWR, _O_RDONLY) != 0)
			throw std::exception("mcScorePhaseSpaceConcentrator:: Cannot open model file to load");
		long buflen = _filelength(ifile);
		void* buf = malloc(buflen);
		if (!buf)
			throw std::exception("can't allocate memory for source model");
		long rlen = _read(ifile, buf, buflen);
		_close(ifile);
		if (rlen != buflen)
			throw std::exception("error in reading source data file");

		source = new mcSourceModelRadialPhsp(srcName.c_str(), nThreads, z0);
		((mcSourceModelRadialPhsp*) source)->readFromMemory(buf);

		free(buf);
	}
	// electron beam from accelerator
	else if (_wcsicmp(radTypeName.c_str(), L"ea_beam") == 0)
	{
		if (fileName.empty())
			throw exception("Please, specify beam particles file name");

		ifstream ibms(XmlParseReaderBase::copyWStringToStlString(fileName.c_str()).c_str());
		if(ibms.fail())
			throw exception("Can not open particles file");
		
		auto src = new mcSourceAcceleratedBeam(srcName.c_str(), nThreads, z0);
		src->loadData(ibms);
		source = src;
	}
	else if (_wcsicmp(radTypeName.c_str(), L"c60") == 0)
	{
		if (!isCoreSetDefined)
			throw exception("Please, indicate radiation source dimensions");
		source = new mcSourceCylindricalC60(srcName.c_str(), nThreads, geomVector3D(x0, y0, z0), geomVector3D(vx, vy, vz), d / 2, h);
	}
	else if (_wcsicmp(radTypeName.c_str(), L"c60_sphere") == 0)
	{
		if (!isCoreSetDefined)
			throw exception("Please, indicate radiation source dimensions");
		source = new mcSourceSphereC60(srcName.c_str(), nThreads, geomVector3D(x0, y0, z0), geomVector3D(vx, vy, vz), radius);
	}
	else if (_wcsicmp(radTypeName.c_str(), L"clinical_electron_beam") == 0)
	{
		source = new mcClinicalElectronBeam(srcName.c_str(), nThreads, energy, z0,
			ewidth, awidth, fsx, fsy);
	}
	else
	{
		if (!isShapeDefined)
			throw exception("Please, indicate radiation source shape");
		
		mc_particle_t particleType = GeometryParser::convert_S2T_ptype(radTypeName.c_str());

		if (_wcsicmp(srcType.c_str(), L"mono") == 0)
		{
			if (_wcsicmp(shapeType.c_str(), L"simple") == 0)
			{
				source = new mcSourceSimpleMono(srcName.c_str(),
					nThreads, particleType, energy, geomVector3D(x0, y0, z0), geomVector3D(vx, vy, vz));
			}
			else if (_wcsicmp(shapeType.c_str(), L"cylindrical") == 0)
			{
			}
			else
				throw exception("Unknown radiation source size");
		}
		else if (_wcsicmp(srcType.c_str(), L"conical") == 0)
		{
			source = new mcSourceAccelerator(srcName.c_str(), nThreads, particleType, energy, z0, size, angle);
		}
		else if (_wcsicmp(srcType.c_str(), L"conical_sigmar") == 0)
		{
			source = new mcSourceXraySigmaR(srcName.c_str(), nThreads, particleType, energy, z0, size, angle);
		}
		else if (_wcsicmp(srcType.c_str(), L"src_leba") == 0)
		{
			source = new mcSourceLEBA(srcName.c_str(), nThreads, sptype, prf_type, energy, ewidth, z0, size);
		}
		else
			throw exception("Unknown radiation source type");
	}

	if (source != nullptr)
	{
		source->setModuleName(sourceModule.c_str());
		source->SetIsStartInside(isStartInside);
	}

	return source;
}
