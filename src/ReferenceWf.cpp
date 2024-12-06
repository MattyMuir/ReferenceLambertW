#include "../include/config.h"
#include "ReferenceWf.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cfenv>

#include <iostream>
#include <format>
#include <numeric>

#define SLEEF_STATIC_LIBS
#include <sleef.h>

#include "rndutil.h"
#include "halley.h"

// (-1/e) rounded towards +Inf
static constexpr float EM_UP = -0.36787942f;

ReferenceWf::ReferenceWf()
{
	mpfr_init2(m, 24);
	mpfr_init2(yLowP0, 24);
	mpfr_init2(yHighP0, 24);
	mpfr_init2(yLowP1, 60);
	mpfr_init2(yHighP1, 60);
}

ReferenceWf::~ReferenceWf()
{
	mpfr_clear(m);
	mpfr_clear(yLowP0);
	mpfr_clear(yHighP0);
	mpfr_clear(yLowP1);
	mpfr_clear(yHighP1);
}

Intervalf ReferenceWf::W0(float x)
{
#if REFERENCEW_STATS
	numEvals++;
#endif

	// Edge cases
	if (x < EM_UP)
		return { NAN, NAN };
	if (x == INFINITY)
		return { FLT_MAX, INFINITY };

	// Save current rounding mode
	int initialRnd = fegetround();

	// === Compute Bracket ===
	// high = ln(x + 1)
	float high;
	if (x > 3.0f)
		high = Sleef_logf_u10(x);
	else
		high = Sleef_log1pf_u10(x);
	high = std::nextafter(high, INFINITY);

	float low;
	if (x > 3)
	{
		// low = ln(x) - ln(ln(x))
		auto [logDown, logUp] = LogUpDown(x);
		logUp = std::nextafter(Sleef_logf_u10(logUp), INFINITY);

		low = sub(logDown, logUp, FE_DOWNWARD);
	}
	else if (x >= 0)
	{
		// low = x / (x + 1)
		float xp1 = add(x, 1, FE_UPWARD);
		low = div(x, xp1, FE_DOWNWARD);
	}
	else
	{
		// low = x * (1 - x * 5)
		float x5 = mul(x, 5, FE_DOWNWARD);
		low = sub(1, x5, FE_UPWARD);
		low = mul(low, x, FE_DOWNWARD);

		// Clamp low above -1
		if (low < -1)
			low = -1;
	}

	// === Halley Iterations ===
	for (size_t i = 0; i < 3; i++)
	{
		low = HalleyW0(x, low, false);
		high = HalleyW0(x, high, true);
	}

	// === Bisection ===
	auto ret = Bisection(x, low, high, true);
	assert(ret.inf == ret.sup || ret.sup == std::nextafterf(ret.inf, INFINITY));

	// Restore rounding mode
	fesetround(initialRnd);

	return ret;
}

Intervalf ReferenceWf::Wm1(float x)
{
#if REFERENCEW_STATS
	numEvals++;
#endif

	// Edge cases
	if (x < EM_UP || x >= 0)
		return { NAN, NAN };

	// Save current rounding mode
	int initialRnd = fegetround();

	// === Compute Bracket ===
	auto [uUp, uDown] = LogUpDown(-x);

	uDown = -uDown;
	uDown = sub(uDown, 1, FE_DOWNWARD);

	uUp = -uUp;
	uUp = sub(uUp, 1, FE_UPWARD);

	// low = -1 - (sqrt(u * 2) + u);
	float low = add(uUp, uUp, FE_UPWARD);
	low = sqrt(low, FE_UPWARD);
	low = add(low, uUp, FE_UPWARD);
	low = sub(-1, low, FE_DOWNWARD);

	// high = -1 - (sqrt(u * 2) + u * 2 / 3)
	float high = add(uDown, uDown, FE_DOWNWARD);
	high = sqrt(high, FE_DOWNWARD);
	uDown = add(uDown, uDown, FE_DOWNWARD);
	uDown = div(uDown, 3, FE_DOWNWARD);
	high = add(high, uDown, FE_DOWNWARD);
	high = sub(-1, high, FE_UPWARD);
	if (!std::isfinite(high))
		high = -1.0f;

	// === Halley Iterations ===
	size_t numIter = ceil(log10(-x) * -0.3 + 1.5);
	for (size_t i = 0; i < 4; i++)
	{
		low = HalleyWm1(x, low, false);
		high = HalleyWm1(x, high, true);
	}

	// === Bisection ===
	auto ret = Bisection(x, low, high, false);
	assert(ret.inf == ret.sup || ret.sup == std::nextafter(ret.inf, INFINITY));

	// Restore rounding mode
	fesetround(initialRnd);

	return ret;
}

#if REFERENCEW_STATS
double ReferenceWf::GetHighPrecRate() const
{
	return (double)numHighPrec / totalBisections;
}

size_t ReferenceWf::GetMaxBisections() const
{
	return maxBisections;
}

double ReferenceWf::GetAvgBisections() const
{
	return (double)totalBisections / numEvals;
}
#endif

Sign ReferenceWf::GetMidpointSign(float x, float midpoint, bool useHighPrec)
{
#if REFERENCEW_STATS
	if (useHighPrec)
		numHighPrec++;
#endif

	mpfr_set_flt(m, midpoint, MPFR_RNDN);

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

	return Sign::Inconclusive;
}

Intervalf ReferenceWf::Bisection(float x, float low, float high, bool increasing)
{
#if REFERENCEW_STATS
	size_t b = 0;
#endif

	for (;;)
	{
#if REFERENCEW_STATS
		b++;
#endif

		// m = (low + high) / 2
		float m = std::midpoint(low, high);

		if (m == low || m == high)
			break; // Bracket cannot be narrowed any further

		// Calculate midpoint sign
		Sign sign = GetMidpointSign(x, m, false);
		if (sign == Sign::Inconclusive)
			sign = GetMidpointSign(x, m, true);

		if (sign == Sign::Inconclusive)
		{
			std::cerr << std::format("Error, ambiguous sign: {}\n", x);
			throw;
		}

		// Update bracket
		if ((sign == Sign::Positive) == increasing)
			high = m;
		else
			low = m;
	}

#if REFERENCEW_STATS
	maxBisections = std::max(maxBisections, b);
	totalBisections += b;
#endif

	return { low, high };
}