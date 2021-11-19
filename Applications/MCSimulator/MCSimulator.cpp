#include "stdafx.h"
#include "GeometryParser.h"
#include "XmlParseReaderBase.h"

#include "../../mc/mc/mcTransport.h"
#include "../../mc/mc/mcScore.h"
#include "../../mc/mc/mcSource.h"
#include "../../mc/mc/mcMedia.h"
#include "../../mc/mc/mcMediumXE.h"
#include "../../mc/mc/mcThread.h"
#include "../../mc/mc/mcVRMLDumper.h"

#include <fstream>
#include <time.h>
#include <ppl.h>
#include <concrtrm.h>

using namespace std;

#define TIME ((double) clock()/CLOCKS_PER_SEC)
#define TIMESINCE(t) ((double) (TIME - (t)))

int _tmain(int argc, _TCHAR* argv [])
{
	if (argc != 3)
	{
		wcout << L"Usage:" << endl;
		wcout << argv[0] << L" params.xml geometry.xml" << endl;
		return -1;
	}

	int i;
	int nThreads = concurrency::GetProcessorCount();
	//int nThreads = 1;

	double simStartTime = TIME;

	mcThread* threads = new mcThread[nThreads];
	for (i = 0; i < nThreads; i++)
		threads[i].setId(i);

	wcout << L"nThreads =  " << nThreads << endl;

	try
	{
		// Используемые среды
		mcMedia media;
		media.addName("H2O700ICRU");
		media.addName("AIR700ICRU");
		media.addName("LUNG700ICRU");
		media.addName("ICRPBONE700ICRU");
		media.addName("CERROBEND700");
		media.addName("PB700ICRU");
		media.addName("W700ICRU");
		media.addName("AL700ICRU");
		media.addName("TI700ICRU");
		media.addName("STEEL700ICRU");
		media.addName("CU700ICRU");
		media.addName("C60");
		media.addName("NAI700ICRU");
		media.addName("SI700ICRU");
		media.addName("MYLAR700ICRU");
		media.addName("HE700ICRU");
		media.addName("170C700ICRU");
		media.addName("226C700ICRU");
		media.addName("BE700ICRU");
		media.addName("TA700ICRU");
		media.addName("AU700ICRU");

		media.initXEFromFile("../data/AcceleratorSimulator.pegs4dat");

		// Pars input files
		XPRNode paramsDoc, geometryDoc;
		XmlParseReaderBase::CreateXPRDocumentFromFile(argv[1], paramsDoc);
		XmlParseReaderBase::CreateXPRDocumentFromFile(argv[2], geometryDoc);

		//
		// Геометрические модули
		//
		mcTransport* tfirst = nullptr;
		mcTransport* tprev = nullptr;

		if (_wcsicmp(geometryDoc.Name.c_str(), L"accelerator") != 0)
			throw exception("Geometry XML file should have root node \"accelerator\"");

		for (auto node : geometryDoc.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"module") == 0)
			{
				mcTransport* t = GeometryParser::ParseTransport(node, &media, nThreads);
				t->setDefaultScore(nThreads);
				if (tprev)
				{
					tprev->setNextTransport(t);
					t->setPreviousTransport(tprev);
				}
				else
					tfirst = t;
				tprev = t;
			}
		}

		//
		// Параметры
		//
		double ecut = 0;
		mcSource* source = nullptr;
		int nHistories = 0;
		int nBanches = 0;
		wstring vrmlFile;
		wstring statFile;
		double trackR = DBL_MAX, trackZ1 = -DBL_MAX, trackZ2 = DBL_MAX, trackEMIN = 0;
		bool doTrack = false, doTrackPhotons = true, doTrackElectrons = true, doTrackPositrons = true;

		if (_wcsicmp(paramsDoc.Name.c_str(), L"input") != 0)
			throw exception("Parameters XML file should have root node \"input\"");

		for (auto node : paramsDoc.Nodes)
		{
			if (_wcsicmp(node.Name.c_str(), L"simulation") == 0)
			{
				for (auto n1 : node.Nodes)
				{
					if (_wcsicmp(n1.Name.c_str(), L"nhistories") == 0)
					{
						nHistories = _wtoi(n1.Text.c_str());
					}
					else if (_wcsicmp(n1.Name.c_str(), L"nbanches") == 0)
					{
						nBanches = _wtoi(n1.Text.c_str());
					}
				}
			}
			//
			// Options
			//
			else if (_wcsicmp(node.Name.c_str(), L"options") == 0)
			{
				for (auto n1 : node.Nodes)
				{
					if (_wcsicmp(n1.Name.c_str(), L"vrmlfile") == 0)
					{
						vrmlFile = n1.Text;
					}
					else if (_wcsicmp(n1.Name.c_str(), L"statfile") == 0)
					{
						statFile = n1.Text;
					}
					else if (_wcsicmp(n1.Name.c_str(), L"trackparticles") == 0)
					{
						doTrack = true;
						for (auto n2 : n1.Nodes)
						{
							if (_wcsicmp(n2.Name.c_str(), L"bounding") == 0)
							{
								for (auto n3 : n2.Nodes)
								{
									if (_wcsicmp(n3.Name.c_str(), L"R") == 0)
										trackR = _wtof(n3.Text.c_str());
									else if (_wcsicmp(n3.Name.c_str(), L"Z1") == 0)
										trackZ1 = _wtof(n3.Text.c_str());
									else if (_wcsicmp(n3.Name.c_str(), L"Z2") == 0)
										trackZ2 = _wtof(n3.Text.c_str());
									else if (_wcsicmp(n3.Name.c_str(), L"EMIN") == 0)
										trackEMIN = _wtof(n3.Text.c_str());
								}
							}
							else if (_wcsicmp(n2.Name.c_str(), L"particles") == 0)
							{
								for (auto n3 : n2.Nodes)
								{
									if (_wcsicmp(n3.Name.c_str(), L"photons") == 0)
										doTrackPhotons = _wcsicmp(n3.Text.c_str(), L"true") == 0 ? true : false;
									else if (_wcsicmp(n3.Name.c_str(), L"electrons") == 0)
										doTrackElectrons = _wcsicmp(n3.Text.c_str(), L"true") == 0 ? true : false;
									else if (_wcsicmp(n3.Name.c_str(), L"positrons") == 0)
										doTrackPositrons = _wcsicmp(n3.Text.c_str(), L"true") == 0 ? true : false;
								}
							}
						}
					}
					else if (_wcsicmp(n1.Name.c_str(), L"transCutoff_elec") == 0)
					{
						for (auto n2 : n1.Nodes)
						{
							if (_wcsicmp(n2.Name.c_str(), L"ecat") == 0)
							{
								ecut = _wtof(n2.Text.c_str());
							}
						}
					}
				}
			}
			//
			// Источник
			//
			else if (_wcsicmp(node.Name.c_str(), L"source") == 0)
			{
				if (source != nullptr)
					throw exception("Only one particles source is supported");
				source = GeometryParser::ParseSource(node, nThreads);
			}
			//
			// Scores
			//
			else if (_wcsicmp(node.Name.c_str(), L"score") == 0)
			{
				mcScore* score = GeometryParser::ParseScore(node, nThreads);
				// Привязка к геометрическому объекту
				mcTransport* t = tfirst;
				while (t)
				{
					if (strcmp(t->getName(), score->getModuleName()) == 0) {
						t->setScore(score);
						break;
					}
					t = (mcTransport*) t->getNextTransport();
				}
				if (t == nullptr)
					throw exception("Geometry module for score not found");
			}
		}

		// Кол-во частиц в одном запуске для каждой нити
		int nParticles = nHistories / nThreads;

		// HACK!!!
		// CUTOFF энергия должна ассоциироваться с модулем.
		// Поскольку там реально не предусмотрена настройка, определяем это глобально здесь.
		// Если соотвествующий параметр не определен или 0, то действуем как прописано в сечениях.
		if (ecut > 0)
		{
			vector<mcMedium*> mm = media.Media();
			for (i = 0; i < (int) mm.size(); i++)
			{
				mcMediumXE* m = (mcMediumXE*) mm[i];
				m->transCutoff_elec = ecut;

				// Test: ручное отключение транспорта электронов в коллиматоре
				//if(m->name_ == "W700ICRU")
				//	m->transCutoff_elec = 100.0;

				// HACK!! До решения вопроса привязки ECUT к модулю принудительно назначаем его
				// для воды и воздуха 0.3 МэВ, как наиболее практичного для транспорта в среде.
				if (m->name_ == "H2O700ICRU" || m->name_ == "AIR700ICRU" ||
					m->name_ == "LUNG700ICRU" || m->name_ == "ICRPBONE700ICRU")
					m->transCutoff_elec = 0.3;
			}
		}

		if (source == nullptr)
			throw exception("Particles ource must be provided");
		bool startinside = source->IsGamma() || source->IsStartInside();
		if (doTrack)
			source->setScoreTrack(trackR, trackZ1, trackZ2, trackEMIN, doTrackPhotons, doTrackElectrons, doTrackPositrons);

		// Запись геометрии в VRML (до симуляции на случай если таковая сломается; 
		// можно бкдет посмотреть простую причину)
		if (!vrmlFile.empty())
		{
			ofstream wfs(vrmlFile.c_str());
			if (wfs.fail())
				throw exception("Can't open wrl file for writing");

			mcVRMLDumper::dumpHead(wfs);
			mcVRMLDumper::dumpWorldAxis(wfs);

			mcTransport* t = tfirst;
			while (t)
			{
				t->dumpVRML(wfs);
				mcScore* s = t->getScore();
				if (s) s->dumpVRML(wfs);

				// Проверяем вложения
				mcTransport* ti = t;
				if (ti != nullptr) ti = ti->getInternalTransport();
				while (ti)
				{
					ti->dumpVRML(wfs);
					mcScore* ss = ti->getScore();
					if (s) ss->dumpVRML(wfs);
					ti = ti->getInternalTransport();
				}

				t = (mcTransport*)t->getNextTransport();
			}
			source->dumpVRML(wfs);

			wcout << L"VRML has been dumped to " << vrmlFile << endl;
		}
		else
			wcout << L"VRML dump skipped" << endl;

		//
		// Симуляция
		//
		mcTransport* tstart = tfirst;
		//if (startinside)
		{
			mcTransport* t = tfirst;
			while (t)
			{
				auto tfound = t->getInternalTransportByName(source->getModuleName());
				if (tfound != nullptr)
				{
					t = tfound;
					break;
				}
				t = (mcTransport*) t->getNextTransport();
			}
			if (t == nullptr)
				throw exception("Gamma source should be related to the object");
			tstart = t;
		}

		double* energySource = new double[nThreads];  // счетчик энергии источника
		for (i = 0; i < nThreads; i++) energySource[i] = 0;

		for (int ib = 0; ib < nBanches; ib++)
		{
			wcout << L"Time:" << TIMESINCE(simStartTime) << L" sec" << L"       Start bunch " << ib + 1 << L" of " << nBanches << endl;

			if (nThreads == 1)
			{
				mcParticle particle;
				for (int ii = 0; ii < nParticles; ii++)
				{
					source->sample(particle, &threads[0]);
					energySource[0] += particle.ke * particle.weight;

					if (startinside)
						tstart->beginTransportInside(particle);
					else
						tstart->beginTransport(particle);
				}
			}
			else
			{
				concurrency::parallel_for(0u, (unsigned)nThreads, [threads, source, tstart, nParticles, energySource, startinside](unsigned it)
				{
					mcParticle particle;
					for (int ii = 0; ii < nParticles; ii++)
					{
						source->sample(particle, &threads[it]);
						energySource[it] += particle.ke * particle.weight;

						if (startinside)
							tstart->beginTransportInside(particle);
						else
							tstart->beginTransport(particle);
					}
				});
			}
		}

		delete [] threads;

		wcout << L"Simulation time =  " << TIMESINCE(simStartTime) << L" sec" << endl;
		wcout << L"Writing results ..." << endl;

		// Сохранение статистики
		double energyDeposited = 0;
		for (i = 0; i < nThreads; i++) energyDeposited += energySource[i];
		wcout << endl;
		wcout << L"Source energy = \t" << energyDeposited << "\tMeV" << endl;

		if (!statFile.empty())
		{
			ofstream sfs(statFile.c_str());
			if (sfs.fail())
				throw exception("Can't open file for writing statistic");

			mcTransport* t = tfirst;
			while (t)
			{
				t->dump(sfs);
				mcScore* s = t->getScore();
				if (s) {
					s->CE2D();
					s->dumpStatistic(sfs);
				}

				// Проверяем вложения
				mcTransport* ti = t;
				if (ti != nullptr) ti = ti->getInternalTransport();
				while (ti)
				{
					ti->dump(sfs);
					mcScore* ss = ti->getScore();
					if (ss) {
						ss->CE2D();
						ss->dumpStatistic(sfs);
					}
					ti = ti->getInternalTransport();
				}

				t = (mcTransport*) t->getNextTransport();
			}

			wcout << L"Statistic has been dumped to " << statFile << endl;
		}
		else
			wcout << L"Statistic dump skipped" << endl;

		// Чистка памяти
		if (source) delete source;

		mcTransport* t = tfirst;
		while (t)
		{
			mcTransport* tt = t;
			t = (mcTransport *) tt->getNextTransport();

			// Проверяем вложения
			mcTransport* ti = tt;
			if (ti != nullptr) ti = ti->getInternalTransport();
			while (ti)
			{
				mcTransport* tii = ti;
				ti = ti->getInternalTransport();
				delete tii;
			}

			delete tt;
		}

		wcout << L"Simulation completed successfuly..." << endl;
	}
	catch (std::exception& e)
	{
		wcout << L"Simulation error: " << e.what() << endl;
	}
	catch (...)
	{
		wcout << L"Unknown simulation error ..." << endl;
	}

	return 0;
}
