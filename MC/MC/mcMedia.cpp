#include "mcMedia.h"
#include "mcMediumXE.h"
#include "mcMediumProton.h"
#include "mcMediumNeutron.h"
#include "mcPhysicsPhoton.h"
#include "mcPhysicsElectron.h"
#include "mcPhysicsPositron.h"
#include "mcPhysicsProton.h"
#include "mcPhysicsNeutron.h"
#include "mcParticle.h"
#include "../geometry/text.h"
#include <fstream>
//#include <string>
//#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

mcMedia::mcMedia(void)
{
	physics_.resize(MCP_NTYPES, nullptr);
	physics_[MCP_PHOTON] = new mcPhysicsPhoton();
	physics_[MCP_NEGATRON] = new mcPhysicsElectron();
	physics_[MCP_POSITRON] = new mcPhysicsPositron();
	physics_[MCP_PROTON] = new mcPhysicsProton();
	physics_[MCP_NEUTRON] = new mcPhysicsNeutron();
}

mcMedia::~mcMedia(void)
{
	int i;
	for (i = 0; i < (int)xes_.size(); i++)
		delete xes_[i];
	for (i = 0; i < (int)protons_.size(); i++)
		delete protons_[i];
	for (i = 0; i < (int)neutrons_.size(); i++)
		delete neutrons_[i];
	for (i = 0; i < (int)physics_.size(); i++)
		delete physics_[i];
}

void mcMedia::addName(const char* mname)
{
	if (!xes_.empty() || !protons_.empty() || !neutrons_.empty())
		throw std::exception("Can't add media names after data initialization");
	mnames_.push_back(mname);
}

short mcMedia::getMediumIdx(const char* mname) const
{
	for (short i = 0; i < (short)mnames_.size(); i++)
		if (mnames_[i] == mname) return i;
	throw std::exception((string("Medium \"") + mname + string("\" does not exist")).c_str());
}

const mcMediumXE* mcMedia::getMediumXE(short idx) const
{
	if (idx >= (short)xes_.size())
		throw std::exception("Medium index for photons an electrons is too big");
	return (mcMediumXE*)xes_[idx];
}

const mcMediumProton* mcMedia::getProtonMedium(short idx) const
{
	if (idx >= (short)protons_.size())
		throw std::exception("Medium index for protons is too big");
	return (mcMediumProton*)protons_[idx];
}

const mcMediumProton* mcMedia::getNeutronMedium(short idx) const
{
	if (idx >= (short)neutrons_.size())
		throw std::exception("Medium index for protons is too big");
	return (mcMediumProton*)neutrons_[idx];
}

void mcMedia::initXEFromStream(istream& is)
{
	if (!xes_.empty())
		throw std::exception("Photon and electron crossections already initialized");
	for (int i = 0; i < (int)mnames_.size(); i++)
		xes_.push_back(new mcMediumXE());

	// „тение данных
	string line, s1, s2;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// ѕровер€ем, нужна ли данна€ среда дл€ загрузки?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				xes_[i]->name_ = line;
				((mcMediumXE*)xes_[i])->read(is);
			}
		}
		getline(is, line, '\n');
	}

	// ѕровер€ем, все ли среды загружены
	string errmedia;
	for (int i = 0; i < (int)xes_.size(); i++)
	{
		if (xes_[i]->status_ != mcMedium::LOADED) {
			errmedia += xes_[i]->name_;
			errmedia += "\n";
		}
	}
	if (!errmedia.empty())
		throw std::exception((string("The following XE media were not loaded succcessfuly:\n") + errmedia).c_str());
}

void mcMedia::initXEFromFile(const string& fname)
{
	ifstream is(fname.c_str());
	if (is.fail())
		throw std::exception((string("Can't open XE data file: ") + fname).c_str());
	initXEFromStream(is);
}

void mcMedia::initProtonDeDxFromStream(istream& is)
{
	if (!protons_.empty())
		throw std::exception("Proton crossectons already initialized");
	int i;
	for (i = 0; i < (int)mnames_.size(); i++)
		protons_.push_back(new mcMediumProton());

	// „тение данных - часть в этой функции полностью аналогична XA, только добавлена проверка версии
	string line, s1, s2, s3, s4;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// ѕровер€ем, нужна ли данна€ среда дл€ загрузки?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				// дополнительно провер€ем версию input file VER=0.0.0
				GetTwoStringsFromLine(s1, s2, s3);
				GetTwoStringsFromLine(s3, s1, s4);
				if ((s2 == "VER") || (s3 == "0.0.0")) {
					protons_[i]->name_ = line;
					((mcMediumProton*)protons_[i])->read(is);
				}
				else {
					//throw std::exception("Wrong Proton media data version"); 
					//в принципе данные могут быть дальше в этом же файле в другой версии, 
					// так что просто не считываем данные
				}
			}
		}
		getline(is, line, '\n');
	}

	// ѕровер€ем, все ли среды загружены
	string errmedia;
	for (int i = 0; i < (int)protons_.size(); i++)
	{
		if (protons_[i]->status_ != mcMedium::LOADED) {
			errmedia += mnames_[i];
			errmedia += "\n";
		}
	}
	if (!errmedia.empty())
		throw std::exception((string("The following Proton media were not loaded succcessfuly:\n") + errmedia).c_str());
}

void mcMedia::initProtonFromFiles(const string& fname, const string& pstardir, const string& icru63dir)
{
	// —тарый вариант тормозных спопосбностей ( остюченко, 2008)
	ifstream is(fname.c_str());
	if (is.fail())
		throw std::exception((string("Can't open Proton data file: ") + fname).c_str());
	initProtonDeDxFromStream(is);

	// TODO: нова€ верси€ тормозных способностей из базы данных PSTAR.
	// Note: эта база содержит весьма ограниченный набор атомов, 
	// по крайней мере в версии с запросами по имени атома л материала.
	// Ќужно разобртьс€, нет ли где-то файлов дл€ всех атомов или это просто формула
	// дл€ которой нужно только брать правильные поенциалы I.
	// 
	// Ќе реализовано! ƒолжно заменить старый вариант.
	//

	// —начала формируем базу данных сечений всех имеющихс€ атомов.
	// ‘айлы дл€ протонов удовлетвор€ют маске "P*.DAT".

	// »дем по всем файлам и формируем базу данных ICRU 63 дл€ протонов.

	//fs::path p = fs::current_path();
	//std::cout << "The current path " << p << " decomposes into:\n"
	//	<< "root-path " << p.root_path() << '\n'
	//	<< "relative path " << p.relative_path() << '\n';

	//std::filesystem::current_path("C:/MCSimulations/protons/");

	/*
	std::vector<FileEntry> CollectFiles(const fs::path & inPath)
	{
		std::vector<fs::path> paths;
		if (fs::exists(inPath) && fs::is_directory(inPath))
		{
			std::filesystem::recursive_directory_iterator dirpos{ inPath };

			std::copy_if(begin(dirpos), end(dirpos), std::back_inserter(paths),
				[](const fs::directory_entry& entry) {
					return entry.is_regular_file();
				}
			);
		}
		std::vector<FileEntry> files(paths.size());
		std::transform(paths.cbegin(), paths.cend(), files.begin(), FileEntry::Create);
		return files;
	}
	*/

	for (const auto& entry : fs::directory_iterator(icru63dir))
	{
		//if (std::regex_match(entry., reg))
		std::cout << entry.path() << std::endl;

	}
}

void mcMedia::initNeutronFromStream(istream& is)
{
	if (!neutrons_.empty())
		throw std::exception("Neutron crossectons already initialized");
	int i;
	for (i = 0; i < (int)mnames_.size(); i++)
		neutrons_.push_back(new mcMediumNeutron());

	// „тение данных - часть в этой функции полностью аналогична XA, только добавлена проверка версии
	string line, s1, s2, s3, s4;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// ѕровер€ем, нужна ли данна€ среда дл€ загрузки?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				// дополнительно провер€ем версию input file VER=0.0.0
				GetTwoStringsFromLine(s1, s2, s3);
				GetTwoStringsFromLine(s3, s1, s4);
				if ((s2 == "VER") || (s3 == "0.0.0")) {
					neutrons_[i]->name_ = line;
					((mcMediumNeutron*)neutrons_[i])->read(is);
				}
				else {
					//throw std::exception("Wrong Neutron media data version"); 
					//в принципе данные могут быть дальше в этом же файле в другой версии, 
					// так что просто не считываем данные
				}
			}
		}
		getline(is, line, '\n');
	}

	// ѕровер€ем, все ли среды загружены
	string errmedia;
	for (int i = 0; i < (int)neutrons_.size(); i++)
	{
		if (neutrons_[i]->status_ != mcMedium::LOADED) {
			errmedia += mnames_[i];
			errmedia += "\n";
		}
	}
	if (!errmedia.empty())
		throw std::exception((string("The following Neutron media were not loaded succcessfuly:\n") + errmedia).c_str());
}

void mcMedia::initNeutronFromFile(const string& fname)
{
	ifstream is(fname.c_str());
	if (is.fail())
		throw std::exception((string("Can't open Neutron data file: ") + fname).c_str());
	initNeutronFromStream(is);
}

const mcPhysics* mcMedia::getPhysics(int ptype) const
{
	if (ptype >= (int)physics_.size())
		throw std::exception("Unsupported particle type");
	return physics_[ptype];
}

const mcMedium* mcMedia::getMedium(int ptype, int idx) const
{
	if (idx >= (int)mnames_.size())
		throw std::exception("Unsupported particle type");
	if (ptype == MCP_PROTON)
		return protons_[idx];
	else if (ptype == MCP_NEUTRON)
		return neutrons_[idx];
	else
		return xes_[idx];
}
