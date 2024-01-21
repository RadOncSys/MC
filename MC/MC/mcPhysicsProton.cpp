#include "mcPhysicsProton.h"
#include "mcMediumProton.h"
#include "mcPhysicsCommon.h"
#include "mcParticle.h"
#include "mcRng.h"
#include "mcThread.h"
#include "mcDefs.h"
#include <float.h>

mcPhysicsProton::mcPhysicsProton(void)
{
}

mcPhysicsProton::~mcPhysicsProton(void)
{
}

bool mcPhysicsProton::Discarge(mcParticle* p, const mcMedium& med, double& edep) const
{
	if (p->ke <= ((const mcMediumProton&)med).transCutoff_proto)
	{
		edep = p->ke;
		DiscardParticle(p);
		return true;
	}
	else
		return false;
}

double mcPhysicsProton::MeanFreePath(double ke, const mcMedium& med, double dens) const
{
	// вместо логарифма используем линейно энергию
	const mcMediumProton& m = (const mcMediumProton&)med;
	double logKE = ke;//log(ke);
	int iLogKE = int(ke);//int(m.iLogKE0_proto + logKE * m.iLogKE1_proto);
	double sigma = (m.sigma0_proto[iLogKE] + logKE * m.sigma1_proto[iLogKE]) * dens;
	return (sigma > 0.0) ? 1 / sigma : DBL_MAX;
}

double mcPhysicsProton::TakeOneStep(mcParticle* p, const mcMedium& med, double& step) const
{
	const mcMediumProton& m = (const mcMediumProton&)med;
	double e_dep = 0;
	// VK
	// Определить размер шага, на который происходит транспорт и занести его в step
	// Определить потери энергии на этом шаге 
	// с учётом страглинга
	// и вернуть их return e_dep 
	// Изменить положение частицы
	// Разыграть рассеяние - изменить направление

	// 1. Определяем средние потери энергии
	// Пока полагаем, что протон на каждом шаге теряет не больше 1% энергии
	e_dep = 0.01 * p->ke;	// Задаём максимальную величину потерь энергии на шаге
	// поправка на убыль энергии вдоль шага
	// берём dE/dx для средней энергии вдоль шага: kE-dE/2.0
	//double mE = kE;		// исходная версия (Отчёт)
	double mE = p->ke - e_dep / 2.0;  // коррекция март 2008
	// вместо логарифма используем линейно энергию
	double logKE = mE; //p->ke;//log(p->ke);
	int iLogKE = int(logKE);//(int) (m.iLogKE0_proto + logKE * m.iLogKE1_proto);
	double dedx = p->regDensityRatio * (m.dedx0_proto[iLogKE] + logKE * m.dedx1_proto[iLogKE]);
	// более короткие шаги (до пересечения с границей или точечного вз. считаем без поправки
	step = MIN(step, e_dep / dedx);
	e_dep = step * dedx;	// если изменился шаг

	// Поправку на то, что пробег больше длины шага не вводим
	// 2. Учитываем страгглинг
	double rel = (p->ke + PMASS) * PMASS_1;
	double sig2 = p->regDensityRatio * m.dEdxStragglingGaussVarianceConstPart_ * step * (1 + SQUARE(rel));

	double r1, r2;
	GaussStandardRnd_by_Marsaglia(p->thread_->rng(), r1, r2); // случайно распределённые по Гауссу величины
	e_dep = e_dep - sqrt(sig2) * r1;
	e_dep = (e_dep < 0) ? 0 : (e_dep > p->ke) ? p->ke : e_dep;

	// 3. Транспортируем частицу
	p->p += p->u * step;
	//p->ke -= e_dep;

	// 4. Разыгрываем новое направление и устанавливаем его
	// Моделируем Мольеровское рассеяние
	double tE = p->ke + PMASS; // полная энергия
	// Не уверен на счёт p->regDensityRatio, но в Nova вроде именно так входит
	double l = p->regDensityRatio * step / m.radLength;
	//double a1	= 1+0.038*log(l); // выражение в скобках
	//double a2	= betasq(PMASS,tE)*tE; // знаменатель b*c*p=E*b^2
	//double a1_2 = a1/a2;
	double a1_2 = (1 + 0.038 * log(l)) / (betasq(PMASS, tE) * tE);
	double th0 = sqrt(184.96 * l * SQUARE(a1_2)); // 184.96=13.6^2	
	// Всё-таки не до конца ясно, а в EGS восстановить не удаётся. Здесь только попытка
	// Моделируем theta и отклонение по
	// см. rpp-2006-book.pdf 27.3 p.262 eq.27.18,27.19
	// а второй угол берём равномерным на 2*pi
	// Надо бы разобраться в том, что в Nova и как его изменить под p+ (масса!)
	// там наверняка лучше, чем здесь
	double theta = fabs(r2 * th0); //eq.27.19
	//смещение игнорируем
	//double sqrtOf12inverse = 0.28867513459481288225457439025098;	
	//double deflection	= r2*th0*t*sqrtOf12inverse+theta*t*0.5;	//eq.27.18
	double phi = p->thread_->rng().rnd() * TWOPI;
	ChangeDirection(cos(theta), sin(theta), cos(phi), sin(phi), p->u);

	// p->mfps = 0; // VK Это зачем? Если мы не разыгрываем каждый раз mfp, а редуцируем его, то эта строка кажестя не верна 

	p->ke -= e_dep;
	return e_dep;
}

double mcPhysicsProton::DoInterruction(mcParticle* p, const mcMedium* med) const
{
	// TODO: Следует серьезнее разобраться с дискретными взаимодействиями.
	// 1. Насколько часто происходит взаимодействие на атомах водорода.
	//    Вероятно их нужно моделировать корректно, так как появляются 2 протона
	//    выделяющие энергию именно в области интереса.
	// 2. В результате ядерных взаимодействий протон меняет направление.
	//    И если не моделировать тяжелые частицы, то хотя-бы определиться с 
	//    энергиями направлениями остаточных протонов, причем в зависимости от их
	//    энергии до взаимодействия.

	// Все дискретные взаимодействия моделируем следующим образом.
	// У протона остается 1/3 кинетической энергии (а куда он движется ???).
	// 1/3 поглощается в точке.
	// 1/3 выносится за пределы области интереса.

	p->ke /= 3.0;
	//cout << ".";

	// Возвращаем энергию, выделившуюся в точке.
	return 2 * p->ke * p->weight;
}
