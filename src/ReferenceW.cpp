#include "../include/ReferenceW.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cfenv>

#include <iomanip>
#include <iostream>
#include <numeric>

#define SLEEF_STATIC_LIBS
#include <sleef.h>

// (-1/e) rounded towards +Inf
static const double EM_UP = -0.3678794411714423;

static inline double add(double x, double y, int rnd)
{
	fesetround(rnd);
	return x + y;
}

static inline double sub(double x, double y, int rnd)
{
	fesetround(rnd);
	return x - y;
}

static inline double mul(double x, double y, int rnd)
{
	fesetround(rnd);
	return x * y;
}

static inline double div(double x, double y, int rnd)
{
	fesetround(rnd);
	return x / y;
}

static inline double sqrt(double x, int rnd)
{
	fesetround(rnd);
	return sqrt(x);
}

static inline std::pair<double, double> ExpUpDown(double x)
{
	double v = Sleef_exp_u10(x);
	return { std::nextafter(v, -INFINITY), std::nextafter(v, INFINITY) };
}

static inline void ExpUpDown(mpfr_t down, mpfr_t up, mpfr_t x)
{
	int isBelow = mpfr_exp(down, x, MPFR_RNDD);
	mpfr_set(up, down, MPFR_RNDN);
	if (isBelow) mpfr_nextabove(up);
}

static inline std::pair<double, double> LogUpDown(double x)
{
	double v = Sleef_log_u10(x);
	return { std::nextafter(v, -INFINITY), std::nextafter(v, INFINITY) };
}

ReferenceW::ReferenceW()
{
	mpfr_init2(m, 53);
	mpfr_init2(yLowP0, 53);
	mpfr_init2(yHighP0, 53);
	mpfr_init2(yLowP1, 150);
	mpfr_init2(yHighP1, 150);
}

ReferenceW::~ReferenceW()
{
	mpfr_clear(m);
	mpfr_clear(yLowP0);
	mpfr_clear(yHighP0);
	mpfr_clear(yLowP1);
	mpfr_clear(yHighP1);
}

Interval ReferenceW::W0(double x)
{
	// Edge cases
	if (x < EM_UP)
		return { NAN, NAN };
	if (x == INFINITY)
		return { DBL_MAX, INFINITY };

	// === Compute Bracket ===
	// high = ln(x + 1)
	double high = Sleef_log1p_u10(x);
	high = std::nextafter(high, INFINITY);

	double low;
	if (x > 3)
	{
		// low = ln(x) - ln(ln(x))
		auto [logDown, logUp] = LogUpDown(x);
		logUp = std::nextafter(Sleef_log_u10(logUp), INFINITY);

		low = sub(logDown, logUp, FE_DOWNWARD);
	}
	else if (x >= 0)
	{
		// low = x / (x + 1)
		double xp1 = add(x, 1, FE_UPWARD);
		low = div(x, xp1, FE_DOWNWARD);
	}
	else
	{
		// low = x * (1 - x * 5)
		double x5 = mul(x, 5, FE_DOWNWARD);
		low = sub(1, x5, FE_UPWARD);
		low = mul(low, x, FE_DOWNWARD);

		// Clamp low above -1
		if (low < -1)
			low = -1;
	}

	// === Halley Iterations ===
	for (size_t i = 0; i < 4; i++)
	{
		low = HalleyW0(x, low, false);
		high = HalleyW0(x, high, true);
	}

	// === Bisection ===
	auto ret = Bisection(x, low, high, true);
	assert(ret.inf == ret.sup || ret.sup == std::nextafter(ret.inf, INFINITY));

	return ret;
}

Interval ReferenceW::Wm1(double x)
{
	// Edge cases
	if (x < EM_UP || x >= 0)
		return { NAN, NAN };

	// === Compute Bracket ===
	auto [uUp, uDown] = LogUpDown(-x);

	uDown = -uDown;
	uDown = sub(uDown, 1, FE_DOWNWARD);

	uUp = -uUp;
	uUp = sub(uUp, 1, FE_UPWARD);

	// low = -1 - (sqrt(u * 2) + u);
	double low = add(uUp, uUp, FE_UPWARD);
	low = sqrt(low, FE_UPWARD);
	low = add(low, uUp, FE_UPWARD);
	low = sub(-1, low, FE_DOWNWARD);

	// high = -1 - (sqrt(u * 2) + u * 2 / 3)
	double high = add(uDown, uDown, FE_DOWNWARD);
	high = sqrt(high, FE_DOWNWARD);
	uDown = add(uDown, uDown, FE_DOWNWARD);
	uDown = div(uDown, 3, FE_DOWNWARD);
	high = add(high, uDown, FE_DOWNWARD);
	high = sub(-1, high, FE_UPWARD);

	// === Halley Iterations ===
	for (size_t i = 0; i < 3; i++)
	{
		low = HalleyWm1(x, low, false);
		high = HalleyWm1(x, high, true);
	}

	// === Bisection ===
	auto ret = Bisection(x, low, high, false);
	assert(ret.inf == ret.sup || ret.sup == std::nextafter(ret.inf, INFINITY));

	return ret;
}

#if TRACK_BISECTIONS
void ReferenceW::LogBisectionStats() const
{
	std::cout << "Num bisections:   " << numBisections << '\n';
	std::cout << "Num inconclusive: " << numInconclusive << '\n';
	std::cout << "Low precision success rate: " << (double)(numBisections - numInconclusive) / numBisections * 100 << "%\n";
}
#endif

ReferenceW::Sign ReferenceW::GetMidpointSign(double x, double midpoint, bool useHighPrec)
{
#if TRACK_BISECTIONS
	if (!useHighPrec) numBisections++;
#endif

	mpfr_set_d(m, midpoint, MPFR_RNDN);

	mpfr_t& yLow = useHighPrec ? yLowP1 : yLowP0;
	mpfr_t& yHigh = useHighPrec ? yHighP1 : yHighP0;

	// Compute exp
	ExpUpDown(yLow, yHigh, m);
	if (mpfr_cmp_ui(m, 0) < 0)
		mpfr_swap(yLow, yHigh);

	// Compute yLow
	mpfr_mul(yLow, yLow, m, MPFR_RNDD);
	mpfr_sub_d(yLow, yLow, x, MPFR_RNDD);

	// Compute yHigh
	mpfr_mul(yHigh, yHigh, m, MPFR_RNDU);
	mpfr_sub_d(yHigh, yHigh, x, MPFR_RNDU);

	int lowCmp = mpfr_cmp_ui(yLow, 0);
	int highCmp = mpfr_cmp_ui(yHigh, 0);

	if (lowCmp >= 0 && highCmp >= 0)
		return Sign::Positive;
	if (lowCmp <= 0 && highCmp <= 0)
		return Sign::Negative;

#if TRACK_BISECTIONS
	if (!useHighPrec) numInconclusive++;
#endif

	return Sign::Inconclusive;
}

Interval ReferenceW::Bisection(double x, double low, double high, bool increasing)
{
	for (;;)
	{
		// m = (low + high) / 2
		double m = std::midpoint(low, high);

		if (m == low || m == high)
			break; // Bracket cannot be narrowed any further

		// Calculate midpoint sign
		Sign sign = GetMidpointSign(x, m, false);
		if (sign == Sign::Inconclusive)
			sign = GetMidpointSign(x, m, true);

		if (sign == Sign::Inconclusive)
		{
			std::cerr << std::setprecision(20);
			std::cerr << "Error, ambiguous sign: " << x << '\n';
			throw;
		}

		// Update bracket
		if ((sign == Sign::Positive) == increasing)
			high = m;
		else
			low = m;
	}

	return { low, high };
}

double ReferenceW::HalleyW0(double x, double w, bool isUpper)
{
	/*
	result			= w - numerator0 / denominator0
	numerator0		= w e^w - x
	denominator0	= e^w (w + 1) - numerator1 / denominator1
	numerator1		= (w + 2)(w e^w - x)
	denominator1	= 2w + 2

	=== Upper Bound Case ===
	- numerator0	is POSITIVE and needs to be rounded DOWN
	- denominator0	is POSITIVE and needs to be rounded UP
	- numerator1	is POSITIVE and needs to be rounded DOWN
	- denominator1	is POSITIVE and needs to be rounded UP

	=== Lower Bound Case ===
	- numerator0	is NEGATIVE and needs to be rounded UP
	- denominator0	is POSITIVE and needs to be rounded UP
	- numerator1	is NEGATIVE and needs to be rounded DOWN
	- denominator1	is POSITIVE and needs to be rounded DOWN
	*/

	int rnd;

	// expUp and expDown
	auto [expDown, expUp] = ExpUpDown(w);

	// wexpDown
	double exp0 = (w > 0) ? expDown : expUp;
	double wexpDown = mul(w, exp0, FE_DOWNWARD);

	// wexpUp
	double exp1 = (w > 0) ? expUp : expDown;
	double wexpUp = mul(w, exp1, FE_UPWARD);

	// numerator0
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	double numerator0 = sub(isUpper ? wexpDown : wexpUp, x, rnd);

	// numerator1
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	double wplus2 = add(w, 2, rnd);
	double numerator1 = sub(wexpDown, x, FE_DOWNWARD);
	numerator1 = mul(numerator1, wplus2, FE_DOWNWARD);

	// denominator1
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	double denominator1 = add(w, w, rnd);
	denominator1 = add(denominator1, 2, rnd);

	// frac1
	double frac1 = div(numerator1, denominator1, FE_DOWNWARD);

	// denominator0
	double denominator0 = add(w, 1, FE_UPWARD);
	denominator0 = mul(denominator0, expUp, FE_UPWARD);
	denominator0 = sub(denominator0, frac1, FE_UPWARD);

	// newW
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	double newW = div(numerator0, denominator0, rnd);
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	newW = sub(w, newW, rnd);

	return newW;
}

double ReferenceW::HalleyWm1(double x, double w, bool isUpper)
{
	/*
	result			= w - numerator0 / denominator0
	numerator0		= w e^w - x
	denominator0	= e^w (w + 1) - numerator1 / denominator1
	numerator1		= (w + 2)(w e^w - x)
	denominator1	= 2w + 2

	=== Upper Bound Case ===
	- numerator0	is NEGATIVE				and needs to be rounded UP
	- denominator0	is NEGATIVE				and needs to be rounded DOWN
	- numerator1	is NEGATIVE when w > -2 and needs to be rounded DOWN
	- denominator1	is NEGATIVE				and needs to be rounded UP when w > -2

	=== Lower Bound Case ===
	- numerator0	is POSITIVE				and needs to be rounded DOWN
	- denominator0	is NEGATIVE				and needs to be rounded DOWN
	- numerator1	is POSITIVE when w > -2 and needs to be rounded DOWN
	- denominator1	is NEGATIVE				and needs to be rounded DOWN when w > -2
	*/

	int rnd;

	// expUp and expDown
	auto [expDown, expUp] = ExpUpDown(w);

	// wexpDown
	double wexpDown = mul(w, expUp, FE_DOWNWARD);

	// wexpUp
	double wexpUp = mul(w, expDown, FE_UPWARD);

	// numerator0
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	double numerator0 = sub(isUpper ? wexpUp : wexpDown, x, rnd);

	// numerator1
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	double wplus2 = add(w, 2, rnd);
	rnd = (wplus2 > 0) ? FE_DOWNWARD : FE_UPWARD;
	double numerator1 = sub((wplus2 > 0) ? wexpDown : wexpUp, x, rnd);
	numerator1 = mul(numerator1, wplus2, FE_DOWNWARD);

	// denominator1
	rnd = (isUpper == (wplus2 > 0)) ? FE_UPWARD : FE_DOWNWARD;
	double denominator1 = add(w, w, rnd);
	denominator1 = add(denominator1, 2, rnd);

	// frac1
	double frac1 = div(numerator1, denominator1, FE_UPWARD);

	// denominator0
	double denominator0 = add(w, 1, FE_DOWNWARD);
	denominator0 = mul(denominator0, expUp, FE_DOWNWARD);
	denominator0 = sub(denominator0, frac1, FE_DOWNWARD);

	// newW
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	double newW = div(numerator0, denominator0, rnd);
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	newW = sub(w, newW, rnd);

	return newW;
}