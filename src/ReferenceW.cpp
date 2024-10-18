#include "ReferenceW.h"

#include <cassert>
#include <cfloat>
#include <cmath>

#include <iomanip>
#include <iostream>

// (-1/e) rounded towards +Inf
static const double EM_UP = -0.3678794411714423;

ReferenceW::ReferenceW()
{
	mpfr_init2(xMpfr, 53);
	mpfr_init2(low, 53);
	mpfr_init2(high, 53);
	mpfr_init2(temp0, 53);
	mpfr_init2(temp1, 53);

	mpfr_init2(m, 53);
	mpfr_init2(yLowP0, 53);
	mpfr_init2(yHighP0, 53);
	mpfr_init2(yLowP1, 150);
	mpfr_init2(yHighP1, 150);
}

ReferenceW::~ReferenceW()
{
	mpfr_clear(xMpfr);
	mpfr_clear(low);
	mpfr_clear(high);
	mpfr_clear(temp0);
	mpfr_clear(temp1);

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

	// Convert x to mpfr
	mpfr_set_d(xMpfr, x, MPFR_RNDN);

	// === Compute Bracket ===
	// high = ln(x + 1)
	mpfr_log1p(high, xMpfr, MPFR_RNDU);

	if (x > 3)
	{
		// low = ln(x) - ln(ln(x))
		mpfr_log(temp0, xMpfr, MPFR_RNDD);

		mpfr_log(temp1, xMpfr, MPFR_RNDU);
		mpfr_log(temp1, temp1, MPFR_RNDU);

		mpfr_sub(low, temp0, temp1, MPFR_RNDD);
	}
	else if (x >= 0)
	{
		// low = x / (x + 1)
		mpfr_add_ui(temp0, xMpfr, 1, MPFR_RNDU);
		mpfr_div(low, xMpfr, temp0, MPFR_RNDD);
	}
	else
	{
		// low = x * (1 - x * 5)
		mpfr_mul_ui(temp0, xMpfr, 5, MPFR_RNDN);
		mpfr_ui_sub(low, 1, temp0, MPFR_RNDU);
		mpfr_mul(low, low, xMpfr, MPFR_RNDD);

		// Clamp low above -1
		if (mpfr_cmp_si(low, -1) < 0)
			mpfr_set_si(low, -1, MPFR_RNDN);
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

	// Convert x to mpfr
	mpfr_set_d(xMpfr, x, MPFR_RNDN);

	// === Compute Bracket ===
	mpfr_neg(temp0, xMpfr, MPFR_RNDN);
	mpfr_log(temp0, temp0, MPFR_RNDU);
	mpfr_neg(temp0, temp0, MPFR_RNDN);
	mpfr_sub_ui(temp0, temp0, 1, MPFR_RNDD);

	mpfr_neg(temp1, xMpfr, MPFR_RNDN);
	mpfr_log(temp1, temp1, MPFR_RNDD);
	mpfr_neg(temp1, temp1, MPFR_RNDN);
	mpfr_sub_ui(temp1, temp1, 1, MPFR_RNDU);

	// low = -1 - (sqrt(u * 2) + u);
	mpfr_mul_2ui(low, temp1, 1, MPFR_RNDU);
	mpfr_sqrt(low, low, MPFR_RNDU);
	mpfr_add(low, low, temp1, MPFR_RNDU);
	mpfr_si_sub(low, -1, low, MPFR_RNDD);

	// high = -1 - (sqrt(u * 2) + u * 2 / 3)
	mpfr_mul_2ui(high, temp0, 1, MPFR_RNDD);
	mpfr_sqrt(high, high, MPFR_RNDD);
	mpfr_mul_2ui(temp0, temp0, 1, MPFR_RNDD);
	mpfr_div_ui(temp0, temp0, 3, MPFR_RNDD);
	mpfr_add(high, high, temp0, MPFR_RNDD);
	mpfr_si_sub(high, -1, high, MPFR_RNDU);

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

ReferenceW::Sign ReferenceW::GetMidpointSign(double x, mpfr_t midpoint, bool useHighPrec)
{
#if TRACK_BISECTIONS
	if (!useHighPrec) numBisections++;
#endif

	mpfr_t& yLow = useHighPrec ? yLowP1 : yLowP0;
	mpfr_t& yHigh = useHighPrec ? yHighP1 : yHighP0;

	// Compute yLow
	mpfr_rnd_t rnd = (mpfr_cmp_ui(midpoint, 0) > 0) ? MPFR_RNDD : MPFR_RNDU;
	mpfr_exp(yLow, midpoint, rnd);
	mpfr_mul(yLow, yLow, midpoint, MPFR_RNDD);
	mpfr_sub_d(yLow, yLow, x, MPFR_RNDD);

	// Compute yHigh
	rnd = (mpfr_cmp_ui(midpoint, 0) > 0) ? MPFR_RNDU : MPFR_RNDD;
	mpfr_exp(yHigh, midpoint, rnd);
	mpfr_mul(yHigh, yHigh, midpoint, MPFR_RNDU);
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

Interval ReferenceW::Bisection(double x, mpfr_t low, mpfr_t high, bool increasing)
{
	for (;;)
	{
		// m = (low + high) / 2
		mpfr_add(m, low, high, MPFR_RNDN);
		mpfr_div_2ui(m, m, 1, MPFR_RNDN);

		if (mpfr_equal_p(m, low) || mpfr_equal_p(m, high))
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
		if (sign == Sign::Positive)
			mpfr_set(increasing ? high : low, m, MPFR_RNDN);
		else
			mpfr_set(increasing ? low : high, m, MPFR_RNDN);
	}

	double lowD = mpfr_get_d(low, MPFR_RNDD);
	double highD = mpfr_get_d(high, MPFR_RNDU);

	return { lowD, highD };
}