#include "rndutil.h"

#include <cfloat>
#include <cmath>
#include <cfenv>

#include <utility>

#define SLEEF_STATIC_LIBS
#include <sleef.h>

float add(float x, float y, int rnd)
{
	fesetround(rnd);
	return x + y;
}

float sub(float x, float y, int rnd)
{
	fesetround(rnd);
	return x - y;
}

float mul(float x, float y, int rnd)
{
	fesetround(rnd);
	return x * y;
}

float div(float x, float y, int rnd)
{
	fesetround(rnd);
	return x / y;
}

float sqrt(float x, int rnd)
{
	fesetround(rnd);
	return sqrtf(x);
}

float fma(float x, float y, float z, int rnd)
{
	fesetround(rnd);
	return fmaf(x, y, z);
}

std::pair<float, float> ExpUpDown(float x)
{
	fesetround(FE_TONEAREST);
	float v = Sleef_expf_u10(x);
	return { std::nextafterf(v, -INFINITY), std::nextafterf(v, INFINITY) };
}

std::pair<float, float> LogUpDown(float x)
{
	fesetround(FE_TONEAREST);
	float v = Sleef_logf_u10(x);
	return { std::nextafterf(v, -INFINITY), std::nextafterf(v, INFINITY) };
}

double add(double x, double y, int rnd)
{
	fesetround(rnd);
	return x + y;
}

double sub(double x, double y, int rnd)
{
	fesetround(rnd);
	return x - y;
}

double mul(double x, double y, int rnd)
{
	fesetround(rnd);
	return x * y;
}

double div(double x, double y, int rnd)
{
	fesetround(rnd);
	return x / y;
}

double sqrt(double x, int rnd)
{
	fesetround(rnd);
	return sqrt(x);
}

float fma(double x, double y, double z, int rnd)
{
	fesetround(rnd);
	return fma(x, y, z);
}

std::pair<double, double> ExpUpDown(double x)
{
	fesetround(FE_TONEAREST);
	double v = Sleef_exp_u10(x);
	return { std::nextafter(v, -INFINITY), std::nextafter(v, INFINITY) };
}

std::pair<double, double> LogUpDown(double x)
{
	fesetround(FE_TONEAREST);
	double v = Sleef_log_u10(x);
	return { std::nextafter(v, -INFINITY), std::nextafter(v, INFINITY) };
}

void ExpUpDown(mpfr_t down, mpfr_t up, mpfr_t x)
{
	int isBelow = mpfr_exp(down, x, MPFR_RNDD);
	mpfr_set(up, down, MPFR_RNDN);
	if (isBelow) mpfr_nextabove(up);
}