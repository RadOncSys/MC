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
#include "mcMendeleev.h"
#include "mcEndfP.h"
#include "mcPStar.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>
#include <ctype.h>

// Флаги вывода служебной информации в процесе отладки
// Вывод осуществляется в корневую папку расчета.
// Имя файла зависит от типа вывода и назначается непосредственно в коде.
// 0 - ничего не выводить
// 1 - вывод сечений реакций протонов
#define MEDIA_TRACE 1

#if MEDIA_TRACE > 0
#include "sstream"
#endif

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

void mcMedia::initProtonFromFiles(const string& pstardir, const string& nuclearDir)
{
#if MEDIA_TRACE > 0
	ofstream fos("ProtonEndfLog.dat");
	fos << "Report about parsing proton database" << endl;
	fos << "----------------------------------------------" << endl;
#endif

	mcMendeleev Table;
	for (int i = 0; i < xes_.size(); i++)
		for (int j = 0; j < xes_[i]->elements_.size(); j++)
			Table.IsNecessary[xes_[i]->elements_[j].atomicNumber - 1] = true;

	// DE/dx from PSTAR
	mcPStar pstardb;
	mcMendeleev ptable(Table);
	pstardb.LoadFromPath(pstardir, ptable);

	// ENDF
	mcEndfDB endfdb;

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
		int Z = 0;

		for (int i = 0; i < AtNum.size(); i++)
		{
			if (isalpha(AtNum[i]))
			{
				AtNum.erase(i, AtNum.size() - i);
				Z = stoi(AtNum);
				break;
			}
		}

		// База данных изотопа для протонов
		if (Table.IsNecessary[Z - 1])
		{
			if(Table.IsLoad[Z - 1])
				throw std::exception((
					string("Something wrong with proton ENDF data! Attempt to load data for element Z= ") +
					to_string(Z) + ", while data alreadt loaded.").c_str());
		
			auto csForElement = std::make_shared<mcEndfNP>();
			csForElement->Load(fs::path(entry.path()).string().c_str(), elementName.c_str());
			csForElement->Z = Z;
			endfdb.Isotopes.push_back(csForElement);
			Table.IsLoad[Z - 1] = true;

#if MEDIA_TRACE == 1
			fos << "Element:\t" << csForElement->ElementName << endl;
			fos << "----------------------------------------------" << endl << endl;
			fos << "Elastic crosssections (MF=3 MT=2)" << endl;
			if (!csForElement->ElasticCrossSections.isEmpty)
				csForElement->ElasticCrossSections.dump(fos);
			else
				fos << endl << "NO Elastic crosssections !!!" << endl << endl;
			fos << endl << endl;
			fos << "Nuclear crosssections (MF=3 MT=5)" << endl;
			if (!csForElement->NuclearCrossSections.isEmpty)
				csForElement->NuclearCrossSections.dump(fos);
			else
				fos << endl << "NO Nuclear crosssections !!!" << endl << endl;
			fos << endl << endl;
#endif
		}
	}

	// Чтобы не мучиться с отладкой в случае проблем сразу проверяем чего не хватает.
	string info;
	for (int i = 0; i < Table.IsLoad.size(); i++)
	{
		if (Table.IsNecessary[i] && !Table.IsLoad[i])
			info += string("Not found proton ENDF element with Z = ") + to_string(i + 1) + "\r\n";
	}
	if(info.size() != 0)
		throw std::exception(info.c_str());

	// Подготавливаем среды для протонов по шаблону EGS
	// и устанавливаем параметры среды в части взаимодействия с электронами (PSTAR)
	if (!protons_.empty())
		throw std::exception("Proton crossectons already initialized");
	for (int i = 0; i < xes_.size(); i++)
	{
		// Конструктор по шаблону EGS устанавливает и копию элементного состава среды
		auto m = new mcMediumProton(*xes_[i]);

		// Тормозные способности по базе данных PSTAR
		m->SetEnergyLoses(pstardb);

		// Упругие и неупругие рассеяния из базы данных ENDF для протонов
		m->SetNuclearCrossSections(endfdb);
		
		m->status_ = mcMedium::LOADED;
		protons_.push_back(m);

#if MEDIA_TRACE > 0
		m->dump(fos);
#endif
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

void mcMedia::initNeutronFromFiles(const string& path, const string& nuclearDir)
{
	ifstream is(path.c_str());
	if (is.fail())
		throw std::exception((string("Can't open Neutron data file: ") + path).c_str());
	initNeutronFromStream(is);

	mcMendeleev Table;

	for (int i = 0; i < xes_.size(); i++)
		for (int j = 0; j < xes_[i]->elements_.size(); j++)
			Table.IsNecessary[xes_[i]->elements_[j].atomicNumber - 1] = true;

	// ENDF
	auto dbData = std::make_shared<std::vector<std::shared_ptr<mcEndfN>>>();

	// Цикл по файлам сечений, в каждом из которых содержатся полные данные для одного изотопа
	for (const auto& entry : fs::directory_iterator(nuclearDir))
	{
		if (!fs::path(entry.path()).has_stem() || !fs::path(entry.path()).has_extension())
			continue;

		//string fname = fs::path(entry.path()).filename().string();
		string fname = fs::path(entry.path()).stem().string();
		string ext = fs::path(entry.path()).extension().string();
		std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

		if (std::toupper(fname[0]) != 'N' || ext != ".DAT")
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
		auto csForElement = std::make_shared<mcEndfN>();
		if (Table.IsNecessary[Z - 1])
		{
			if (Table.IsLoad[Z - 1])
				throw std::exception((
					string("Something wrong with neutron ENDF data! Attempt to load data for element Z= ") +
					to_string(Z) + ", while data alreadt loaded.").c_str());
			csForElement->Load(fs::path(entry.path()).string().c_str(), elementName.c_str());
			dbData->push_back(csForElement);
			Table.IsLoad[Z - 1] = true;
		}
	}

	string info;
	for (int i = 0; i < Table.IsLoad.size(); i++)
	{
		if (Table.IsNecessary[i] && !Table.IsLoad[i])
			info += string("Not found neutron ENDF element with Z = ") + to_string(i + 1) + "\r\n";
	}
	if (info.size() != 0)
		throw std::exception(info.c_str());

	initNeutronCSFromVector(dbData);
}

void mcMedia::initNeutronCSFromVector(std::shared_ptr<std::vector<std::shared_ptr<mcEndfN>>> dbData)
{
	for (int i = 0; i < neutrons_.size(); i++)
	{
		((mcMediumNeutron*)neutrons_[i])->ENDFdata = dbData;
		((mcMediumNeutron*)neutrons_[i])->createNDB();
	}
}

const mcPhysics* mcMedia::getPhysics(int ptype) const
{
	if (ptype >= (int)physics_.size())
		throw std::exception("Unsupported particle type");
	return physics_[ptype];
}

const mcMedium* mcMedia::getMedium(int ptype, int idx) const
{
	if (ptype == MCP_PHOTON || ptype == MCP_NEGATRON || ptype == MCP_POSITRON)
	{
		if (idx >= (int)xes_.size())
			throw std::exception("EGS media index out of range");
		return xes_[idx];
	}
	else if (ptype == MCP_PROTON)
	{
		if (idx >= (int)protons_.size())
			throw std::exception("Proton media index out of range");
		return protons_[idx];
	}
	else if (ptype == MCP_NEUTRON)
	{
		if (idx >= (int)neutrons_.size())
			throw std::exception("Neutron media index out of range");
		return neutrons_[idx];
	}
	else
		throw std::exception("Unsupported particle type");
}
