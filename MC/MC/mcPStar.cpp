#include "mcPStar.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

mcPStarRecord mcPStarTable::GetDataForEnergy(double e)
{
    mcPStarRecord record;
    if (Data.size() > 1)
    {
        for (int i = 1; i < Data.size(); i++)
        {
            if (Data[i].D[0] >= e)
            {
                record.D[0] = e;
                double f = (log(e) - log(Data[i - 1].D[0])) / (log(Data[i].D[0]) - log(Data[i - 1].D[0]));
                record.D[1] = Data[i - 1].D[1] * (1 - f) + Data[i].D[1] * f;
                record.D[2] = Data[i - 1].D[2] * (1 - f) + Data[i].D[2] * f;
                record.D[3] = Data[i - 1].D[3] * (1 - f) + Data[i].D[3] * f;
                break;
            }
        }
    }
    return record;
}

void mcPStarTable::LoadFromStream(std::istream& is)
{
    // Пустые строки просто пропускаем
    // Начало таблицы после строки с подписями колонок

    std::getline(is, InternaName, '\n');
    TrimLine(InternaName);

    std::string line;
    std::getline(is, line, '\n');
    bool isTableFound = false;
    mcPStarRecord record;

    while (!is.fail())
    {
        std::vector<std::string> ss;
        if (line.size() != 0 && TrimLine(line) != 0 && GetStringArray(line, ss, "\t ") == 7)
        {
            if (!isTableFound && ss[6] == "DETOUR")
            {
                isTableFound = true;
            }
            else if(isTableFound)
            {
                record.D[0] = atof(ss[0].c_str());
                record.D[1] = atof(ss[1].c_str());
                record.D[2] = atof(ss[2].c_str());
                record.D[3] = atof(ss[3].c_str());
                Data.push_back(record);
            }
        }
        std::getline(is, line, '\n');
    }
}

mcPStar::mcPStar()
{
}

mcPStar::~mcPStar()
{
}

std::shared_ptr<mcPStarTable> mcPStar::GetTableForElement(int Z) const
{
    for (int i = 0; i < Tables.size(); i++)
    {
        if (Tables[i]->Z == Z)
            return Tables[i];
    }
    return nullptr;
}

void mcPStar::LoadFromPath(const std::string& path, mcMendeleev& table)
{
    for (int i = 0; i < table.IsNecessary.size(); i++)
    {
        if (!table.IsNecessary[i]) continue;
        
        bool isTbaleFound = false;

        for (const auto& entry : fs::directory_iterator(path))
        {
            if (!fs::path(entry.path()).has_stem() || !fs::path(entry.path()).has_extension())
                continue;

            std::string ext = fs::path(entry.path()).extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
            if (ext != ".DAT")
                continue;

            std::vector<std::string> fna;
            std::string fname = fs::path(entry.path()).stem().string();
            if (GetStringArray(fname, fna, "_") != 3)
                continue;
            int Z = atoi(fna[1].c_str());
            if (Z != (i + 1)) continue;

            std::ifstream is(entry.path());
            if (is.fail())
                throw std::exception((std::string("Can't open PStar table file: ") + entry.path().string()).c_str());

            auto table = std::make_shared<mcPStarTable>();
            table->ElementName = fna[0];
            table->Z = Z;
            table->A = atoi(fna[2].c_str()); 

            table->LoadFromStream(is);

            Tables.push_back(table);
            isTbaleFound = true;
        }

        if (!isTbaleFound)
        {
            char str[256];
            sprintf_s(str, 256, "PSTAR table for Z = %i not found", (i + 1));
            throw std::exception((const char*)str);
        }
    }
}
