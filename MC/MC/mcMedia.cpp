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
//#include "mcCSNuclear.h"
#include "mcEndfP.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>
#include <ctype.h>

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

	// Чтение данных
	string line, s1, s2;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// Проверяем, нужна ли данная среда для загрузки?
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

	// Проверяем, все ли среды загружены
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

	// Чтение данных - часть в этой функции полностью аналогична XA, только добавлена проверка версии
	string line, s1, s2, s3, s4;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// Проверяем, нужна ли данная среда для загрузки?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				// дополнительно проверяем версию input file VER=0.0.0
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

	// Проверяем, все ли среды загружены
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

void mcMedia::initProtonFromFiles(const string& fname, const string& nuclearDir)
{
	// Старый вариант тормозных спопосбностей (Костюченко, 2008)
	ifstream is(fname.c_str());
	if (is.fail())
		throw std::exception((string("Can't open Proton data file: ") + fname).c_str());
	initProtonDeDxFromStream(is);
	

	//TablicaMendeleeva H = Loaded O = Loaded
	//for (media1:media_last)
	//protons_->elements_->
	Mendeleev Table;
	Table.init();

	for (int i = 0; i < xes_.size(); i++)
		for (int j = 0; j < xes_[i]->elements_.size(); j++)
		{
			if (Table.isNecessary[xes_[i]->elements_[j].atomicNumber] == false)
				Table.isNecessary[xes_[i]->elements_[j].atomicNumber] = true;
		}



	// Объект, в который сначала закачиваем всю баз данных сечений
	//auto dbData = std::make_unique<std::vector<std::unique_ptr<mcCSNuclear>>>();

	// ICRU-63
	//auto dbData = std::make_unique<std::vector<mcCSNuclear>>();

	// ENDF
	vector<mcEndfP> dbData;

	// Цикл по файлам сечений, в каждом из которых содержатся полные данные для одного изотопа
	for (const auto& entry : fs::directory_iterator(nuclearDir))
	{
		if (!fs::path(entry.path()).has_stem() || !fs::path(entry.path()).has_extension())
			continue;

		//string fname = fs::path(entry.path()).filename().string();
		string fname = fs::path(entry.path()).stem().string();
		string ext = fs::path(entry.path()).extension().string();
		std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

		if (std::toupper(fname[0]) != 'P' || ext != ".DAT")
			continue;

		// Метку атомного элемента берем из имени файла.
		string elementName = std::string(&fname[2]);
		string AtNum = elementName;
		int Z;

		for (int i = 0; i < AtNum.size(); i++)
		{
			if (isalpha(AtNum[i]))
			{
				AtNum.erase(i, AtNum.size() - i);
				Z = stoi(AtNum);
				break;
			}
		}

		// База данных изотопа
		//mcCSNuclear csForElement;
		mcEndfP csForElement;
		if (Table.isNecessary[Z])
		{
			csForElement.Load(fs::path(entry.path()).string().c_str(), elementName.c_str());
			dbData.push_back(csForElement);
		}
		//ifstream isIcru(entry.path().c_str());
	}
	initProtonCSFromVector(&dbData);
}

void mcMedia::initProtonCSFromVector(std::vector<mcEndfP>* dbData)
{
	for (int i = 0; i < protons_.size(); i++)
	{
		((mcMediumProton*)protons_[i])->ENDFdata = *dbData;
		((mcMediumProton*)protons_[i])->createDB();
	}
}

void mcMedia::initNeutronFromStream(istream& is)
{
	if (!neutrons_.empty())
		throw std::exception("Neutron crossectons already initialized");
	int i;
	for (i = 0; i < (int)mnames_.size(); i++)
		neutrons_.push_back(new mcMediumNeutron());

	// Чтение данных - часть в этой функции полностью аналогична XA, только добавлена проверка версии
	string line, s1, s2, s3, s4;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// Проверяем, нужна ли данная среда для загрузки?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				// дополнительно проверяем версию input file VER=0.0.0
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

	// Проверяем, все ли среды загружены
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

void Mendeleev::init()
{
	for (int i = 0; i < 119; i++)
	{
		isNecessary.push_back(false);
		isLoad.push_back(false);
	}
}
