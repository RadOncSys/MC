#include "mcSamplers.h"


double mcSamplers::SamplePrism(double rnd)
{
	return 2.0 * (rnd - 0.5);
}

double mcSamplers::SampleTriangle(double rnd)
{
	return (rnd < 0.5) ? 2.0 * (sqrt(0.5 * rnd) - 0.5) : 2.0 * (0.5 - sqrt(0.5 * (1.0 - rnd)));
}

double mcSamplers::SampleGauss(double rnd)
{
	// Комулятивным распределением последнего является функция ошибок.
	// Самплится обратная функция ошибка, апроксимация которой описана на:
	// http://en.wikipedia.org/wiki/Error_function

	const double a = 0.140012;
	const double pa = 2.0 / (PI * a);

	rnd = 2.0 * (rnd - 0.5);
	double sign = rnd < 0 ? -1.0 : 1.0;
	rnd *= rnd;

	double f1 = log(1.0 - rnd);
	double f2 = pa + 0.5*f1;

	// Между распределением гаусса и функцией ошибок.
	// Гаусс ~exp(-x^2/2)
	// Erf ~exp(-x^2).
	// Мы работаем Гауссом.
	double f = sqrt(2.0 *(sqrt(f2*f2 - f1 / a) - f2));

	return sign * f;
}

double mcSamplers::SampleGauss2D(double rnd)
{
	const double s = sqrt(2.0);
	return s * ((rnd < 0.5) ? -sqrt(-log(rnd * 2.0)) : sqrt(-log((1.0 - rnd) * 2.0)));
}

std::pair<double, double> mcSamplers::SampleGaussBoxMuller(double rnd1, double rnd2)
{
	const double factor = sqrt(-2 * log(rnd1));
	double val1 = factor * cos(2 * PI * rnd2);
	double val2 = factor * sin(2 * PI * rnd2);

	if (rnd1 < 0.5) 
	{
		val2 = -val2;
		
	}

	if (rnd2 < 0.5)
	{
		val1 = -val1;
	}

	return { val1, val2 };
}

double mcSamplers::SampleExponent2D(double C)
{
	double sign = 1.0;
	if (C < 0.5)
	{
		C *= 2.0;
	}
	else
	{
		sign = -1.0;
		C = 2.0 * (1.0 - C);
	}

	// Соотвествует очень большому расстоянию.
	// Исключается из расчета, так как может создать проблемы.
	if (C < MINDELTA)
		return 1.0 / MINDELTA;

	// На каждом шаге решаем задачу пересечени параболы с прямой линией

	// Стартовая точка
	double y0 = 1;
	double x0 = log(y0 / C);
	double xstep = 1;

	while (xstep > 0.001)
	{
		xstep = x0;
		double A = C * exp(x0);
		double a = A * 0.5;
		double b = A - 1.0;
		double c = b - x0;
		x0 += (sqrt(b*b - 4.0*a*c) - b) / A;
		y0 = 1.0 + x0;
		x0 = log(y0 / C);
		xstep = fabs(x0 - xstep);
	}

	return sign * x0;
}
