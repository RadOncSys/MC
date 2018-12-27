#include "mcScorePhaseSpaceConcentrator.h"
#include "mcParticle.h"
#include "mcThread.h"
#include "mcTransport.h"
#include "mcSourceModelRadialPhsp.h"

#include <vector>
using namespace std;

mcScorePhaseSpaceConcentrator::	mcScorePhaseSpaceConcentrator(
	const char* module_name, int nThreads,
	enum mc_particle_t ptype, bool isXray, double focus,
	int ne, int nr, int nthet, int naxial, 
	double emax, double rmax, double thetmax, 
	const char* stat_file)
:mcScore(module_name, nThreads)
,ptype_(ptype)
,ne_(ne)
,nr_(nr), nthet_(nthet)
,naxial_(naxial), focus_(focus)
,emax_(emax)
,rmax_(rmax)
,model_file_(stat_file)
{
	isxray_ = isXray ? 1 : 0;
	de_ = emax_ / ne_;
	dr_ = 2.0 / nr_;

	// Координатой поворота здесь является не угол, а величина 1 - cos(theta)
	// В этом случае мы сможем избавить от синусов и косинусов при самплинге.
	// Частицы в обратном направлении не рассматриваются
	thetmax_abs_ = thetmax * PI / 180.0;
	thetmax_ = 1 - cos(thetmax_abs_);
	if(thetmax_ >= 1.0)
		throw std::exception("mcScorePhaseSpaceConcentrator:: возможность движения частиц в обратном направлении не рассматривается");
	dthet_ = thetmax_ / nthet_;
	daxial_ = PI / naxial_;

	nnt_ = nthet_ * naxial_;
	nnr_ = nr_ * nnt_;
	nne_ = (isXray ? ne_ : ne_ + 2) * nnr_;
	data_ = new double[nne_ * nThreads_];
	memset(data_, 0, nne_ * nThreads_ * sizeof(double));
}

mcScorePhaseSpaceConcentrator::~mcScorePhaseSpaceConcentrator()
{
	if(data_ != nullptr)
		delete [] data_;
}

void mcScorePhaseSpaceConcentrator::ScoreFluence(const mcParticle& particle)
{
	if(particle.t != ptype_)
		return;

	int iThread = particle.thread_->id();
	etotal_[iThread] += particle.ke * particle.weight;

	int idx = particleIdx(particle);
	if(idx >= 0)
		data_[iThread * nne_ + idx] += particle.weight;
}

int mcScorePhaseSpaceConcentrator::particleIdx(const mcParticle& particle) const
{
	// Переносим частицу в плоскость самплинга распределения
	geomVector3D p = particle.p + (particle.u * (-particle.p.z()/particle.u.z()));

	int eidx;
	eidx = int(particle.ke / de_);

	if(isxray_ == 0)
	{
		if(particle.ke == 1.17) 
			eidx = ne_;
		else if(particle.ke == 1.33) 
			eidx = ne_ + 1;
		else if(eidx >= ne_) 
			return -1;
	}
	else if(eidx >= ne_) 
		return -1;

	double r = p.lengthXY();
	// !! Вероятно при каких то почти горизонтальных движениях возможно переполнение
	// integer, при котором индекс становится равным 0 и частицы попадают 
	// в нулевой бин по радиусу, если их не отсечь дополнительно по радиусу.
	if(r >= rmax_) return -1;
	r = 1.0  - cos(PI * r / rmax_);
	int ridx = int(r / dr_);
	if(ridx < 0 || ridx >= nr_) return -1;

	// Векторы ортогональной системы координат с осями
	// vp - вдоль линии, соединяющея фокус с положением частицы
	geomVector3D vp(p.x(), p.y(), focus_);
	vp.normalize();

	//
	// Отклонение от фокуса
	//
	// GG 20150112 Изменение модели для угла отклонения от направления из фокуса.
	// Угловое распределение тормозных пучков и рассеянных даже для C60 имеет 
	// ярко выраженный экспоненциальный характер даже не смотря 
	// на квадратичную зависимость телесного угла от планарного. 
	// Как реализовать экспоненциональное преобразование координат не ясно
	// из-за неопределенности показателя экспоненты.
	// Поэтому просто используем линейное масштабирование на весь диапазон
	// в надежде, что при очень малых углах экспонента то же линейна и 
	// большая часть фотонов именно при малых углах.

	double thet;
	if (eidx >= ne_)
	{
		// Нерассеяные фотоны регистрируются с разрешением в 4 раза выше
		thet = (1 - (vp * particle.u)) * 4.0;
		if (thet >= thetmax_)
			return -1;
		thet = sin(PI * 0.5 * thet / thetmax_);
	}
	else
	{
		thet = acos(vp * particle.u) / thetmax_abs_;
		if (thet >= 1)
			return -1;
	}

	// Не понятно почему, но по аналогии с радиусом есть проблема усиления нулевого бина, 
	// хотя и менее значительная. Следующая строка ее убирает.
	if(thet < 0)
		return -1;
	int thetidx = int(thet * nthet_);
	if(thetidx >= nthet_) 
		return -1;

	// Азимут
	geomVector3D vu = vp - particle.u;
	vu(2) = 0;
	vu.normalize();
	vp(2) = 0;
	vp.normalize();
	double a = acos(vp * vu);
	int aidx = int(a / daxial_);
	if(aidx >= naxial_) 
		return -1;

	return eidx * nnr_ +  ridx * nnt_ + thetidx * naxial_ + aidx;
}

void mcScorePhaseSpaceConcentrator::dumpStatistic(ostream& os) const
{
	mcScore::dumpStatistic(os);
	
	os << "ptype_:\t" << ptype_ << endl;
	os << "isxray_:\t" << isxray_ << endl;
	os << "focus_:\t" << focus_ << endl;
	os << "ne_:\t" << ne_ << endl;
	os << "nr_:\t" << nr_ << endl;
	os << "nthet_:\t" << nthet_ << endl;
	os << "naxial_:\t" << naxial_ << endl;
	os << "emax_:\t" << emax_ << endl;
	os << "rmax_:\t" << rmax_ << endl;
	os << "thetmax_abs_:\t" << thetmax_abs_ << endl;
	os << "thetmax_:\t" << thetmax_ << endl;
	os << "model_file_:\t" << model_file_.c_str() << endl;

	// Суммируем данные потоков
	double i, count = 0;
	double* dtmp = new double[nne_];
	memset(dtmp, 0, nne_ * sizeof(double));
	double* psrc = data_;
	double* pdest = dtmp;
	for(i=0; i< nne_; i++, pdest++, psrc++)
	{
		for(int it=0; it<nThreads_; it++)
			*pdest += psrc[it * nne_];
		count += *pdest;
	}
	os << "count = \t" << count << endl;

	// Генерируем файл модели для самплинга частиц из данной модели фазового пространства
	mcSourceModelRadialPhsp modelFile("Histogram based source", 1, 0.0);
	modelFile.createFromDistribution(isxray_, ne_, nr_, nthet_, naxial_, 
		emax_, rmax_, thetmax_abs_ * 180.0 / PI, focus_, dtmp);

	int len;
	void* memobj = modelFile.saveToMemory(len);
	if(memobj == nullptr)
		throw std::exception("mcScorePhaseSpaceConcentrator:: Cannot create model object in memory");

	FILE* modelfile;
	if(fopen_s(&modelfile, model_file_.c_str(), "wb") != 0)
		throw std::exception("mcScorePhaseSpaceConcentrator:: Cannot open model file for writing");
	fwrite(memobj, 1, len, modelfile);
	fclose(modelfile);

	// Для исключения проблем на этапе сохранения и восстановления сохраняем модель и восстанавливаем
	mcSourceModelRadialPhsp model("Histogram based source", 1, 0.0);
	model.readFromMemory(memobj);

	free(memobj);

	// Самплим частицы из модели
	mcRng rng;
	rng.init(33, 37);
	mcParticle particle;
	double* dmodel = new double[nne_];
	memset(dmodel, 0, nne_ * sizeof(double));

	for(i=0; i<100000000; i++)
	{
		model.sample(particle, rng);
		int idx = particleIdx(particle);
		if(idx >= 0)
			dmodel[idx] += particle.weight;
	}

	// Интегральные характеристики потоков
	dumpDistributions(os, dtmp, dmodel);

	os << endl << "M O D E L" << endl << endl;
	os << model;

	delete [] dtmp;
	delete [] dmodel;
}

// Нормировка распределения следующим образом:
// каждый элемент массива делится на сумму всех элементов и умножается на 10^6.
void normalizeDistribution(std::vector<double>& d)
{
	double f = 0;
	int i;
	for(i=0; i<(int)d.size(); i++) f += d[i];
	f = 1e6 / f;
	for(i=0; i<(int)d.size(); i++) d[i] *= f;
}

void mcScorePhaseSpaceConcentrator::dumpDistributions(ostream& os, double* dsim, double* dmodel) const
{
	// Интегральные характеристики потоков
	int ne = (isxray_ ? ne_ : ne_ + 2);
	vector<double> espec(ne, 0), espec_2(ne, 0);
	vector<double> rspec(nr_, 0), rspec_2(nr_, 0);
	vector<double> tspec(nthet_, 0), tspec_2(nthet_, 0);
	vector<double> aspec(naxial_, 0), aspec_2(naxial_, 0);

	vector<double> r1spec(nr_, 0), r1spec_2(nr_, 0);
	vector<double> t1spec(nthet_, 0), t1spec_2(nthet_, 0);
	vector<double> a1spec(naxial_, 0), a1spec_2(naxial_, 0);

	vector<double> r2spec(nr_, 0), r2spec_2(nr_, 0);
	vector<double> t2spec(nthet_, 0), t2spec_2(nthet_, 0);
	vector<double> a2spec(naxial_, 0), a2spec_2(naxial_, 0);

	for(int ie=0; ie < ne_; ie++)
	{
		double se = 0, se_2 = 0;
		for(int ir=0; ir < nr_; ir++)
		{
			double st = 0, st_2 = 0;
			for(int it=0; it < nthet_; it++)
			{
				double sa = 0, sa_2 = 0;
				for(int ia=0; ia < naxial_; ia++)
				{
					double f = dsim[ie * nnr_ +  ir * nnt_ + it * naxial_ + ia];
					aspec[ia] += f;
					sa += f;
					f = dmodel[ie * nnr_ +  ir * nnt_ + it * naxial_ + ia];
					aspec_2[ia] += f;
					sa_2 += f;
				}
				tspec[it] += sa;
				st += sa;
				tspec_2[it] += sa_2;
				st_2 += sa_2;
			}
			rspec[ir] += st;
			se += st;
			rspec_2[ir] += st_2;
			se_2 += st_2;
		}
		espec[ie] = se;
		espec_2[ie] = se_2;
	}

	// Две линии кобальта если применимо
	if(isxray_ == 0)
	{
		for(int ir=0; ir < nr_; ir++)
		{
			for(int it=0; it < nthet_; it++)
			{
				for(int ia=0; ia < naxial_; ia++)
				{
					double f = dsim[ne_ * nnr_ +  ir * nnt_ + it * naxial_ + ia];
					espec[ne_] += f;
					r1spec[ir] += f;
					t1spec[it] += f;
					a1spec[ia] += f;
					f = dsim[(ne_+1) * nnr_ +  ir * nnt_ + it * naxial_ + ia];
					espec[ne_+1] += f;
					r2spec[ir] += f;
					t2spec[it] += f;
					a2spec[ia] += f;

					f = dmodel[ne_ * nnr_ +  ir * nnt_ + it * naxial_ + ia];
					espec_2[ne_] += f;
					r1spec_2[ir] += f;
					t1spec_2[it] += f;
					a1spec_2[ia] += f;
					f = dmodel[(ne_+1) * nnr_ +  ir * nnt_ + it * naxial_ + ia];
					espec_2[ne_+1] += f;
					r2spec_2[ir] += f;
					t2spec_2[it] += f;
					a2spec_2[ia] += f;
				}
			}
		}
	}

	// Нормировка всех распределений
	normalizeDistribution(espec);
	normalizeDistribution(espec_2);
	normalizeDistribution(rspec);
	normalizeDistribution(rspec_2);
	normalizeDistribution(tspec);
	normalizeDistribution(tspec_2);
	normalizeDistribution(aspec);
	normalizeDistribution(aspec_2);
	if (isxray_ == 0)
	{
		normalizeDistribution(r1spec);
		normalizeDistribution(r1spec_2);
		normalizeDistribution(t1spec);
		normalizeDistribution(t1spec_2);
		normalizeDistribution(a1spec);
		normalizeDistribution(a1spec_2);
		normalizeDistribution(r2spec);
		normalizeDistribution(r2spec_2);
		normalizeDistribution(t2spec);
		normalizeDistribution(t2spec_2);
		normalizeDistribution(a2spec);
		normalizeDistribution(a2spec_2);
	}

	int i;
	os << endl << "Распределение по энергии (МэВ)" << endl;
	os << "--------------------------------------" << endl;
	for(i=0; i<ne; i++)
		os << (i+0.5) * de_ << "\t" << espec[i] << "\t" << espec_2[i] << endl;

	os << endl << "Распределение по радиусу (см)" << endl;
	os << "-------------------------------------" << endl;
	for (i = 0; i < nr_; i++)
	{
		os << rmax_ * acos(1 - dr_ * (i + 0.5)) / PI << "\t" << rspec[i];
		if (isxray_ == 0)
			os << "\t" << r1spec[i] << "\t" << r2spec[i];
		os << "\t" << rspec_2[i];
		if (isxray_ == 0)
			os << "\t" << r1spec_2[i] << "\t" << r2spec_2[i];
		os << endl;
	}

	os << endl << "Распределение по углу (градусы)" << endl;
	os << "---------------------------------------" << endl;
	for (i = 0; i < nthet_; i++)
	{
		os << (i + 0.5) * (thetmax_abs_ / nthet_) * 180.0 / PI << "\t" << tspec[i];
		os << "\t" << tspec_2[i];
		os << endl;
	}

	if (isxray_ == 0)
	{
		os << endl << "Распределение по углу (нерассеянные фотоны гамма-источника)" << endl;
		os << "---------------------------------------------------------------" << endl;
		for (i = 0; i < nthet_; i++)
		{
			os << acos(1.0 - asin((i + 0.5) * dthet_)) * 180.0 / PI << "\t" << tspec[i];
			os << "\t" << t1spec[i] << "\t" << t2spec[i] << "\t" << t1spec_2[i] << "\t" << t2spec_2[i];
			os << endl;
		}
	}

	os << endl << "Распределение по азимуту (градусы)" << endl;
	os << "------------------------------------------" << endl;
	for (i = 0; i < naxial_; i++)
	{
		os << ((i + 0.5) * daxial_) * 180.0 / PI << "\t" << aspec[i];
		if (isxray_ == 0)
			os << "\t" << a1spec[i] << "\t" << a2spec[i];
		os << "\t" << aspec_2[i];
		if (isxray_ == 0)
			os << "\t" << a1spec_2[i] << "\t" << a2spec_2[i];
		os << endl;
	}

}

void mcScorePhaseSpaceConcentrator::dumpVRML(ostream& os) const
{
	os << "# mcScorePhaseSpaceConcentrator Score: " << name_ << endl;
	if(transport_ == nullptr)
	{
		os << "# Transport not set. Dump not possible!" << endl;
		return;
	}
	const geomMatrix3D& mttow = transport_->MT2W();

	int ir, it, count = 0;
	int da = 15; // шаг по углу 15 градусов
	double mPi = PI / 180;

	os << "Shape {" << endl;
	os << "  appearance Appearance {" << endl;
	os << "    material Material {" << endl;
	os << "      emissiveColor " << red_ << ' ' << green_ << ' ' << blue_ << endl;
	os << "    }" << endl;
	os << "  }" << endl;
	os << "  geometry IndexedLineSet {" << endl;

	os << "    coord Coordinate {" << endl;
	os << "      point [" << endl;

	// Концентрические круги
	for(ir=1; ir<=nr_; ir++) {
		double r = (rmax_ / PI) * acos(1 - dr_ * ir);
		for(it=0; it<360; it+=da) {
			geomVector3D p = geomVector3D(r*sin(mPi*it), r*cos(mPi*it), 0) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			p = geomVector3D(r*sin(mPi*(it+da)), r*cos(mPi*(it+da)), 0) * mttow;
			os << "        " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
			count++;
		}
	}

	os << "      ]" << endl;
	os << "    }" << endl;

	os << "    coordIndex [" << endl;
	for(it=0; it<count; it++)
		os << "      " << 2*it << ' ' << 2*it+1 << " -1" << endl;
	os << "    ]" << endl;
	os << "  }" << endl;
	os << "}" << endl;
}
