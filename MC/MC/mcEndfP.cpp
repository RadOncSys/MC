#include "mcEndfP.h"
#include "../geometry/text.h"
#include <fstream>
#include <filesystem>

using namespace std;

double mcEndfRecord::ParseValue(const char* s, int n)
{
	double f = atof(s);
	int i = 0;
	for (; i < n; i++)
		if (s[i] == '+' || s[i] == '-') break;
	if (i < n)
		f *= pow(10.0, atof(s + i));
	return f;
}

void mcEndfP::Load(const char* fname, const char* ename)
{
	// Separators
	static const string beginSeparator("*** C O N T E N T S ***");
	
	ifstream isEndf(fname);
	if (isEndf.fail())
		throw exception((string("Can't open Proton data file: ") + ename).c_str());

	// Читаем строки текста одну за другой и выбираем нужную информацию
	string line, s1, s2, s3, s4;
	std::getline(isEndf, line, '\n');

	// Состояния указыват в каком месте парсинга мы находимся и потому как интерпитируем строки
	bool isInData = false;
	int pointCount = 0;

	mcEndfRecord record;

	while (!isEndf.fail())
	{
		if(line.size() < 80)
			throw exception((string("Wrong ENDF line length ") + line).c_str());
		memcpy(&record, line.c_str(), 80);

		// Начало новой энергии
		if (!isInData)
		{
			if (line.find(beginSeparator) != string::npos)
				isInData = true;
		}

		// Последняя строка файла. Прерываем не дожидаясь ошибки.
		else if (record.Stblt[0] == '-' && record.Stblt[1] == '1')
			break;

		// Длина строки из девяток является разделителем между таблицами
		else if (string(record.LineNumber, 5) == "99999")
		{
			pointCount = 0;
		}

		// Используем только MF=3 (сечения реакций) / MT=5 (сумма всех реакций за исключением отдельно оговоренных)
		// и     MF=6 (энерго-угловые распределени) / MT=5

		// TODO: Проверить предыдущую гипотезу по другим атомам.
		// Т.е. на предмет того что в других изотопах нет отдельно (от MT=5) описанных реакций.

		// Сечения
		else if (record.MF[0] == ' ' && record.MF[1] == '3' && 
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')
		{
			if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '2')
			{
				CrossSections.ninterpolations = atoi(record.c[4]);
				if (CrossSections.ninterpolations > 1)
					throw exception("UUnexpected multiple interpolation types");
			}
			else if (record.LineNumber[3] == ' ' && record.LineNumber[4] == '3')
			{
				CrossSections.npoints = atoi(record.c[0]);
				CrossSections.interpolationType = atoi(record.c[1]);
				CrossSections.Energies.resize(CrossSections.npoints, 0);
				CrossSections.Values.resize(CrossSections.npoints, 0);
			}
			else
			{
				for (int ii = 0; ii < 6; ii+=2)
				{
					if (pointCount < CrossSections.npoints)
					{
						CrossSections.Energies[pointCount] = mcEndfRecord::ParseValue(record.c[ii], 11);
						CrossSections.Values[pointCount] = mcEndfRecord::ParseValue(record.c[ii+1], 11);
						pointCount++;
					}
				}
			}
		}

		// Энерго-угловые распределения
		else if (record.MF[0] == ' ' && record.MF[1] == '6' && 
			record.MT[0] == ' ' && record.MT[1] == ' ' && record.MT[2] == '5')
		{
		}

		getline(isEndf, line, '\n');
	}
}

void mcEndfP::Clear()
{
	//Energies.clear();
}
