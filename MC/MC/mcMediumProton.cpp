#include "mcMediumProton.h"
#include "../geometry/text.h"
#include <iostream>
#include "mcPhysicsCommon.h"
#include "mcEndfP.h"
//#include <math.h>
//���������� ������ ���������� ������� ��������������� � ��������, ��������:
#ifndef InverseRadiationLength
#define InverseRadiationLength InverseRadiationLength_DahlApproximation
#endif InverseRadiationLength

// �������� �������� ������������ ����� �������� 
// c ������� ������ A [�/����],
// � ������� ���� (������� �������) Z (� �������� ������ ���������)
// ����������� � ����������� Dahl'a
// rpp-2006-book.pdf 27.4.1 p.264 (eq.27.22)
// ��������� � �������� Tsai'� (27.20) � ��������� ����� 2.5%, �� ����������� ����� (5%)
// � 1/(g/cm^2)
double InverseRadiationLength_DahlApproximation(const double A, const double Z)
{
	return Z * (Z + 1) * log(287 / sqrt(Z)) / (716.4 * A);
}

// �������� �������� ������������ ����� �������� ���������� �� n ���������.
// ��� i-�� �������� 0<=i<n
// w[i] - ������������� ��������(?) ��� 
// A[i] - ������� ����� A [�/����],
// Z[i]	- ����� ���� (������� �����)
// rpp-2006-book.pdf 27.4.1 p.263 (eq.27.23)
// ��� ������� ������������ ����� ���������� �������� ��������
// InverseRadiationLength (�������� �����������)
// � 1/(g/cm^2)
double InverseRadiationLength(const double* A, const double* Z, const double* w, const int n)
{
	double s = 0;
	for (int i = 0; i < n; i++)
		s += InverseRadiationLength(A[i], Z[i]) * w[i];
	return s;
}

string CLEARFROMALPHA(string x);

//Необходимо задать конкретную формулу соответствующим с макросом, например:
#ifndef InverseRadiationLength
#define InverseRadiationLength InverseRadiationLength_DahlApproximation
#endif InverseRadiationLength

// Величина обратная радиационной длине вещества 
// c атомной массой A [г/моль],
// и зарядом ядра (атомным номером) Z (в единицах заряда электрона)
// вычисленная в приближении Dahl'a
// rpp-2006-book.pdf 27.4.1 p.264 (eq.27.22)
// совпадает с формулой Tsai'я (27.20) с точностью лучше 2.5%, за исключением гелия (5%)
// в 1/(g/cm^2)
double InverseRadiationLength_DahlApproximation(const double A, const double Z)
{
	return Z * (Z + 1) * log(287 / sqrt(Z)) / (716.4 * A);
}

// Величина обратная радиационной длине вещества состоящего из n элементов.
// Для i-го элемента 0<=i<n
// w[i] - относительный массовый(?) вес 
// A[i] - атомная масса A [г/моль],
// Z[i]	- заряд ядра (атомный номер)
// rpp-2006-book.pdf 27.4.1 p.263 (eq.27.23)
// Для расчёта радиационной длины отдельного элемента вызывает
// InverseRadiationLength (принятое приближение)
// в 1/(g/cm^2)
double InverseRadiationLength(const double* A, const double* Z, const double* w, const int n)
{
	double s = 0;
	for (int i = 0; i < n; i++)
		s += InverseRadiationLength(A[i], Z[i]) * w[i];
	return s;
}

// Возвращает rms радиус ядра в фм (1E-13 см)
// см. Wilson, et al 1991 rp1257.pdf 4.5.2 (eq.4.84,4.85) с исправлениями:
// предполагаем, что в формуле 4.84 ошибка - лишний корень
// "приведённое выражение, предполагает, что распределение ядерной материи, 
// является функцией Гаусса. Такое предположение годится для лёгких ядер (light-weight) 
// и менее употребимо  для At >> 20"
// [VK add 06.08]
double	rmsNuclearRadius(int At)
{
	// ядерный формфактор
	double ac =
		(At == 1) ? 0.84 :
		(At == 2) ? 2.71 :
		(At == 3) ? 1.78 :
		(At == 4) ? 1.63 :
		((At >= 6) && (At <= 14)) ? 2.4 :
		(0.82 * pow(At, 1.0 / 3.0) + 0.58);	// At>=16

	// см также rpp-2006-book.pdf p.71
	// протон Charge Radius = 0.875+-0.007 фм
	// Правда не ясно, почему через A, а не через Z у Уильямса - с нейтроном не сойдётся.

	// Предполагаем, что в формуле 4.84 ошибка - лишний корень
	return sqrt(SQUARE(ac) - 0.64);
}

// ИСПОЛЬЗУЕМАЯ реально функция в 1-ой законченной версии
// Полное сечение ядерных взаимодействий по модели Трипати для лёгких систем
// Результат см2.
// Ap, Zp - атомный номер и заряд налетающей частицы (ядра)
// At, Zt - атомный номер и заряд ядра-мишени
// KE - кинетическая энергия налетающей частицы (ядра) в МэВ
double sigmaTripathiLight(int Ap, int Zp, int At, int Zt, double KE) 
{
	double r0 = 100.0 * 1.1E-15; // в сантиметрах, чтобы сразу возвращать рез. в см2
	double Z3 = 1.0 / 3.0;
	double Ap3 = pow(Ap, Z3);
	double At3 = pow(At, Z3);
	double Apt3 = Ap3 + At3;
	double Ecm = ECenterOfMass(At * PMASS, Ap * PMASS, Ap * PMASS + KE); // приближённое вычисление масс
	double Ecm3 = pow(Ecm, Z3);
	double SL = 1.2 + 1.6 * (1 - exp(-KE / 15.0));
	double X1 = 2.83 - 3.1E-2 * At + 1.7E-4 * SQUARE(At); // для протонов
	double Xm = 1.0 - X1 * exp(-KE / (X1 * SL));
	double D = 1.85 + 0.16 / (1.0 + exp((500.0 - KE) / 200.0));
	double T1 = 18; // для протонов, хотя есть подозрение, что 18 
	double CE = D * (1.0 - exp(-KE / T1)) - 0.292 * exp(-KE / 792.0) * cos(0.229 * pow(KE, 0.453));// было ошибочно 0.291,
	double S = Ap3 * At3 / Apt3;
	double deltaE = (1.85 + 0.16 / Ecm3) * S - CE + 0.91 * (At - 2 * Zt) * Zp / (At * Ap);
	//double rTrms	= At3;//1.3*At3; // root mean square radius in fm (по аналогии с ф. Шена)
	//double rPrms	= Ap3;//1.3*Ap3; // моя аппроксимация значений и предположение о единицах
	double rTrms = rmsNuclearRadius(At); // root mean square radius in fm 
	double rPrms = rmsNuclearRadius(Ap); // аппроксимация Willson et al. 1991 [VK add 06.08]
	double rT = 1.29 * rTrms;
	double rP = 1.29 * rPrms;
	double R = rP + rT + 1.2 * (Apt3) / Ecm3;
	double B = 1.44 * Zp * Zt / R;
	double RC =
		(At == 1 && Zt == 1) ? 7.0 : // p+p - моё предположение - линейная экстарполяция. Если это всё вообще применимо для pp
		(At == 2 && Zt == 1) ? 13.5 : // p+d
		(At == 3 && Zt == 2) ? 21.0 : // p+3He
		(At == 4 && Zt == 2) ? 27.0 : // p+4He
		(At == 6 && Zt == 3) ? 2.2 : // p+4He
		1.0; // по аналогии с другими формулами где нет этого множителя и общей тенденцией
	double fff = r0 * (Apt3 + deltaE);
	double sigmaTL = PI * SQUARE(fff) * (1 - RC * B / Ecm) * Xm;
	sigmaTL = (sigmaTL > 0.0) ? sigmaTL : 0.0;
	return sigmaTL;
}

double mcMediumProton::microsigmaforelement(int A, int Z, double kE) const
{
	double SIGMA = 0.0;
	kE *= 1000000;
	bool isFound = false;
	int i = 0;
	string elName = to_string(Z);
	if (A < 10)
		elName += "00" + to_string(A);
	else if (A < 100)
		elName += "0" + to_string(A);
	else elName += to_string(A);
	for (i = 0; i < ENDFdata.size(); i++)
	{
		if (CLEARFROMALPHA(ENDFdata[i].ElementName) == elName)
		{
			isFound = true;
			break;
		}
	}
	if (!isFound)
		return SIGMA;	//���� ������ �� ������ � ���� ������ ENDF ������������ 0
	//throw exception((string("Nucleus with ID: ") + elName + string(" was not found.")).c_str());
	if (ENDFdata[i].NuclearCrossSections.isEmpty)
		return SIGMA;	//���� ��� ������ �� MF=3 MT=5 ������������ 0
	if (kE <= ENDFdata[i].NuclearCrossSections.Energies[0])
		return SIGMA;
	SIGMA = ENDFdata[i].NuclearCrossSections.get_sigma(kE);
	return SIGMA / pow(10,24);
}

double sigmaENDF(int A, int Z, int kE, vector<mcEndfP>* ENDF)
{
	double SIGMA = 0.0;
	kE *= 1000000;
	bool isFound = false;
	int i = 0;
	string elName = to_string(Z);
	if (A < 10)
		elName += "00" + to_string(A);
	else if (A < 100)
		elName += "0" + to_string(A);
	else elName += to_string(A);
	for (i = 0; i < ENDF->size(); i++)
	{
		if (CLEARFROMALPHA(ENDF->at(i).ElementName) == elName)
		{
			isFound = true;
			break;
		}
	}
	if (!isFound)
		return SIGMA;	//���� ������ �� ������ � ���� ������ ENDF ������������ 0
		//throw exception((string("Nucleus with ID: ") + elName + string(" was not found.")).c_str());
	if (ENDF->at(i).NuclearCrossSections.isEmpty)
		return SIGMA;	//���� ��� ������ �� MF=3 MT=5 ������������ 0
	if (kE <= ENDF->at(i).NuclearCrossSections.Energies[0])
		return SIGMA;
	else
	{
		for (int j = 0; j < ENDF->at(i).NuclearCrossSections.Energies.size(); j++)
			if (kE < ENDF->at(i).NuclearCrossSections.Energies[j])
			{
				if (j == 0)
					break;
				SIGMA = ENDF->at(i).NuclearCrossSections.Values[j - 1] + 
					(kE - ENDF->at(i).NuclearCrossSections.Energies[j - 1]) *
					(ENDF->at(i).NuclearCrossSections.Values[j] - ENDF->at(i).NuclearCrossSections.Values[j - 1]) /
					(ENDF->at(i).NuclearCrossSections.Energies[j] - ENDF->at(i).NuclearCrossSections.Energies[j - 1]);
				break;
			}
	}
	return SIGMA;
}

mcMediumProton::mcMediumProton(void)
	:transCutoff_proto(1.0)
{
}

mcMediumProton::~mcMediumProton(void)
{
}

// Атомный вес среды, г/моль SUM{Ni*Ai}
const double mcMediumProton::AtomicWeight() const
{
	double A = 0.0; // Атомный вес
	for (vector<mcElement>::const_iterator el = elements_.begin(); el != elements_.end(); el++) 
	{
		A += el->atomicMass * el->partsByNumber;
	};
	return A;
}

//--------------------------------
// Генерация данных (физика!)
//--------------------------------

//----------------------------------------------------------------------------------
// Параметры страгглинга (разброса) dE/dx на некотором пути в Гауссовом приближении
//----------------------------------------------------------------------------------

// генерирует постоянную (по энергии и пути) часть вариации Гауссова приближения разброса dE/dx. 
// т.е. {sigma^2 of dE/dx} / [path *  (1-(totalE/Mass)^2)]
// исходная формула взяты из диссертации Н.М.Соболевского
const double mcMediumProton::gdEdxStragglingGaussVarianceConstPart()
{
	//(sigma^2 of dE/dx) 
	// = const				* f1(path)	* f2(projectile)		* f3(element)
	// = {0.3*(EMASS)^2}	* {path}	* {1+(totalE/Mass)^2}	* {density*w*Z*A}
	// => for medium from more tham 1 elements (VK) =>
	// = {0.3*(EMASS)^2}	* {path}	* {1+(totalE/Mass)^2}	* {density*SUMi(w[i]Z[i]/A[i])
	// = {this function return}	* {path} * {1+(totalE/Mass)^2}
	// {this function return} = const * f3(medium) = {0.3*(EMASS)^2} * {density*SUMi(w[i]Z[i]/A[i])

	// учитываем, что Wi=Ni/A
	// зависимая от состава среды часть
	dEdxStragglingGaussVarianceConstPart_ = 0.0;
	for (vector<mcElement>::iterator el = elements_.begin(); el != elements_.end(); el++) {
		dEdxStragglingGaussVarianceConstPart_ += el->partsByNumber * el->atomicNumber; // wi*Zi/Ai=ni*Zi/A
	}

	dEdxStragglingGaussVarianceConstPart_ *= 0.3 * SQUARE(EMASS) * density_ / AtomicWeight(); //0.3*SQUARE(EMASS)/A=0.07833601179626508/A
	return dEdxStragglingGaussVarianceConstPart_;
};

// генерирует и возвращает величину радиационной длины сложного вещества 
// rpp-2006-book.pdf 27.4.1 p.263 (eq.27.23)
// для расчёта радиационной длины отдельного элемента вызывает InverseRadiationLength 
// (принятое приближение (в версии 2007 года это приближение Dahl'а))
// в см!!!
const double mcMediumProton::gRadiationLength()
{
	radLength = 0.0;
	for (vector<mcElement>::iterator el = elements_.begin(); el != elements_.end(); el++) {
		radLength += InverseRadiationLength(el->atomicMass, el->atomicNumber) * el->partsByNumber * el->atomicMass;
	}
	radLength = AtomicWeight() / (radLength * density_);
	return radLength;
}


// Вычисляет коэффициенты линейной аппроксимации для каждого диапазона s по двум точкам ax+b,
void coeff_calc(const vector<double>& s, vector<double>& a, vector<double>& b)
{
	int n = (int)s.size();
	a.assign(n, 0.0); b.assign(n, 0.0); // очищаем вектора, устанавливаем размер
	for (int i = 0; i < n - 1; i++) {
		a[i] = s[i + 1] - s[i]; // делим на (i+1) - (i)
		b[i] = s[i];
	}
	//экстраполяция в область больших значений - с коэффициентами предыдущей ячейки
	a[n - 1] = a[n - 2];
	b[n - 1] = b[n - 2];
	return;
}

// mfp=1/(S) 
// даёт sigma*dens*Na/A [1/cm]
// для налетающей частицы с массой (в единицах массы протона) Ap 
// и зарядом (в единицах заряда электрона) Zp
// По умолчанию для протона (Ap = Zp = 1)
const void mcMediumProton::gSigmaInelastic(int Ap, int Zp)
{
	double S;
	vector<double>sigma_in;

	for (int i = 0; i < kEmax(); i++) {
		S = 0.0; // длина свободного пробега
		for (vector<mcElement>::iterator el = elements_.begin(); el != elements_.end(); el++) {
			S += sigmaTripathiLight(Ap, Zp,
				ROUND(el->atomicMass), ROUND(el->atomicNumber),
				i + 1) * el->partsByNumber;
		}
		//mfp_in_1_[i]=S*density_*NAVOGADRO/AtomicWeight();
		sigma_in.push_back(S * NAVOGADRO * density_ / AtomicWeight()); // mfp=1/(sigma_in)
	}
	// Не оптимизмруем, чтобы не запутаться, вычисляем коэффициенты во втором проходе
	coeff_calc(sigma_in, sigma1_proto, sigma0_proto);
}

void mcMediumProton::read(istream& is)
{
	const double distanceUnit = 1.0;
	status_ = FAILED;

	string line, s1, s2;
	vector<double> a;
	vector<int> ia;

	// Line
	getline(is, line, '\n');
	if (is.fail()) return;

	// Тип среды, плотность
	string eletype;
	GetTwoStringsFromLine(line, eletype, s2);
	line = s2;
	this->density_ = atof(ParseLine(line, "RHO").c_str());

	// Элементарный состав
	int i, ne = atoi(ParseLine(line, "NE").c_str());
	this->elements_.resize(ne);
	for (i = 0; i < ne; i++)
	{
		// Line
		getline(is, line, '\n');
		if (is.fail()) return;

		mcElement& e = this->elements_[i];
		strcpy_s(e.atomicSymbol, 3, ParseLine(line, "ASYM").c_str());
		e.atomicNumber = atoi(ParseLine(line, "Z").c_str());
		e.atomicMass = atof(ParseLine(line, "A").c_str());
		e.partsByNumber = atof(ParseLine(line, "PZ").c_str());
	}

	// Выше стандарт из Nova
	// Больше нам пока ничего не надо, кроме таблицы, из которой надо взять dE/dx
	// Это уже таблица нашего типа
	for (;;) { // Ищем начало таблицы TABLE
		// Line
		getline(is, line, '\n');
		if (is.fail()) return;
		if (line.substr(0, 5) == "TABLE")break;
	}
	// 1 строка - имена переменных
	// 2 строка - размерности
	// 3 строка и далее - данные, до END
	// По идее надо иметь возможность опускать 1 и вторую строки
	getline(is, line, '\n');
	if (is.fail()) return;
	vector<string> ss; int nl;
	nl = GetStringArray(line, ss, "\t"); // разбираем заголовок на столбцы
	int ikE = -1, idEdx = -1;
	for (int i = 0; i < nl; i++) { // Ищем столбец с энергией и столбец с dEdx
		if (ss[i] == "kE") { ikE = i; }
		if (ss[i] == "dE/dx") { idEdx = i; }
	}
	// В общем случае надо просто переходить к следующей таблице. 
	// Сейчас полагаем, что таблица одна и если данных в ней нет, то всё
	if ((ikE == -1) || (idEdx == -1)) { throw std::exception("kE or dE/dx data not found"); }
	getline(is, line, '\n'); // строка с размерностями - не проверяем, просто пропукскаем
	vector<double>kE, dEdx;
	getline(is, line, '\n');
	i = 0;
	while (line.substr(0, 3) != "END") {// проверяем конец таблицы или нет
		vector<double> data;
		int nd = GetFloatArray(line, data);
		if (nd != nl) { throw std::exception("Wrong Table"); }
		// В файле данных dE/dx в МэВ*см^2/г, адалее будем работать в линейных единицах, 
		// поэтому сразу перводим, домножая на плотность, которую уже считали
		kE.push_back(data[ikE]); dEdx.push_back(data[idEdx] * density_); //забираем данные
		getline(is, line, '\n');
		i++;
	}
	// Сейчас предполагаем, что в файле данных kE=[1,2,...?], в противном случае данные отвергаем
	for (i = 0; i < kE.size(); i++) {
		if (kE[i] != double(i + 1)) { throw std::exception("kE is not [1,2,...]"); }
	};
	coeff_calc(dEdx, dedx1_proto, dedx0_proto); // вычисляем коэффициенты для линейной интерполяции

	// Данные загрузили, но надо ещё и расчитать недостающие
	gdEdxStragglingGaussVarianceConstPart();
	gRadiationLength();
	gSigmaInelastic();
	gRadiationLength();

	status_ = LOADED;
}

void mcMediumProton::createDB()
{
	double S;
	vector<double>sigma_endf;
	vector<double>sigma_;

	for (int i = 0; i < kEmax(); i++) {
		S = 0.0; // длина свободного пробега
		for (vector<mcElement>::iterator el = elements_.begin(); el != elements_.end(); el++) {
			S += sigmaENDF(ROUND(el->atomicMass), ROUND(el->atomicNumber), i, &ENDFdata)/pow(10,24) * el->partsByNumber;
		}
		//mfp_in_1_[i]=S*density_*NAVOGADRO/AtomicWeight();
		sigma_endf.push_back(S * NAVOGADRO * density_ / AtomicWeight()); // mfp=1/(sigma_in)
		sigma_.push_back(S); // mfp=1/(sigma_in)
	}
	// Не оптимизмруем, чтобы не запутаться, вычисляем коэффициенты во втором проходе
	coeff_calc(sigma_endf, sigma1_proto, sigma0_proto);
}

string CLEARFROMALPHA (string x)
{
	for (int i = 0; i < x.length(); i++)
		if (x[i] > '9')
			x.erase(i, 1);
	return x;
}