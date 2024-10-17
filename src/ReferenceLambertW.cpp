#include "ReferenceLambertW.h"

#include <cassert>
#include <iostream>
#include <iomanip>

#include <mpfr.h>

// (-1/e) rounded towards +Inf
static const double EM_UP = -0.3678794411714423;

Interval Bisection(double x, mpfr_t low, mpfr_t high, bool increasing, mpfr_prec_t precision)
{
	mpfr_t m;
	mpfr_init2(m, 53);

	mpfr_t yLow, yHigh;
	mpfr_init2(yLow, precision);
	mpfr_init2(yHigh, precision);
	for (;;)
	{
		// m = (low + high) / 2
		mpfr_add(m, low, high, MPFR_RNDN);
		mpfr_div_2ui(m, m, 1, MPFR_RNDN);

		if (mpfr_equal_p(m, low) || mpfr_equal_p(m, high))
			break; // Bracket cannot be narrowed any further

		// Compute yLow
		mpfr_rnd_t rnd = (mpfr_cmp_ui(m, 0) > 0) ? MPFR_RNDD : MPFR_RNDU;
		mpfr_exp(yLow, m, rnd);
		mpfr_mul(yLow, yLow, m, MPFR_RNDD);
		mpfr_sub_d(yLow, yLow, x, MPFR_RNDD);

		// Compute yHigh
		rnd = (mpfr_cmp_ui(m, 0) > 0) ? MPFR_RNDU : MPFR_RNDD;
		mpfr_exp(yHigh, m, rnd);
		mpfr_mul(yHigh, yHigh, m, MPFR_RNDU);
		mpfr_sub_d(yHigh, yHigh, x, MPFR_RNDU);

		int lowCmp = mpfr_cmp_ui(yLow, 0);
		int highCmp = mpfr_cmp_ui(yHigh, 0);

		if (lowCmp >= 0 && highCmp >= 0)
			mpfr_set(increasing ? high : low, m, MPFR_RNDN);
		if (lowCmp <= 0 && highCmp <= 0)
			mpfr_set(increasing ? low : high, m, MPFR_RNDN);

		if (lowCmp < 0 && highCmp > 0)
		{
			std::cerr << std::setprecision(20);
			std::cerr << "Error, ambiguous sign: " << x << '\n';
			throw;
		}
	}

	double lowD = mpfr_get_d(low, MPFR_RNDD);
	double highD = mpfr_get_d(high, MPFR_RNDU);

	mpfr_clear(m);
	mpfr_clear(yLow);
	mpfr_clear(yHigh);

	return { lowD, highD };
}

Interval ReferenceW0(double x)
{
	static const mpfr_prec_t W0Precision = 107;

	// Edge cases
	if (x < EM_UP)
		return { NAN, NAN };
	if (x == INFINITY)
		return { DBL_MAX, INFINITY };

	// Initialize variables
	mpfr_t xMpfr, low, high;
	mpfr_init2(xMpfr, 53);
	mpfr_init2(low, 53);
	mpfr_init2(high, 53);

	// Convert x to mpfr
	mpfr_set_d(xMpfr, x, MPFR_RNDN);

	// === Compute Bracket ===
	// high = ln(x + 1)
	mpfr_log1p(high, xMpfr, MPFR_RNDU);

	if (x > 3)
	{
		// low = ln(x) - ln(ln(x))
		mpfr_t logXDown;
		mpfr_init2(logXDown, 53);
		mpfr_log(logXDown, xMpfr, MPFR_RNDD);

		mpfr_t logXUp;
		mpfr_init2(logXUp, 53);
		mpfr_log(logXUp, xMpfr, MPFR_RNDU);
		mpfr_log(logXUp, logXUp, MPFR_RNDU);

		mpfr_sub(low, logXDown, logXUp, MPFR_RNDD);

		mpfr_clear(logXDown);
		mpfr_clear(logXUp);
	}
	else if (x >= 0)
	{
		// low = x / (x + 1)
		mpfr_t xp1;
		mpfr_init2(xp1, 53);
		mpfr_add_ui(xp1, xMpfr, 1, MPFR_RNDU);
		mpfr_div(low, xMpfr, xp1, MPFR_RNDD);

		mpfr_clear(xp1);
	}
	else
	{
		// low = x * (1 - x * 5)
		mpfr_t x5;
		mpfr_init2(x5, 53);
		mpfr_mul_ui(x5, xMpfr, 5, MPFR_RNDN);
		mpfr_ui_sub(low, 1, x5, MPFR_RNDU);
		mpfr_mul(low, low, xMpfr, MPFR_RNDD);

		// Clamp low above -1
		if (mpfr_cmp_si(low, -1) < 0)
			mpfr_set_si(low, -1, MPFR_RNDN);

		mpfr_clear(x5);
	}

	// === Bisection ===
	auto ret = Bisection(x, low, high, true, W0Precision);
	assert(ret.sup == std::nextafter(ret.inf, INFINITY));

	mpfr_clear(xMpfr);
	mpfr_clear(low);
	mpfr_clear(high);

	return ret;
}

// Initial lower/upper bounds are from
// Chatzigeorgiou, "Bounds on the Lambert Function and Their Application to the Outage Analysis of User Cooperation"
Interval ReferenceWm1(double x)
{
	static const mpfr_prec_t Wm1Precision = 82;

	// Edge cases
	if (x < EM_UP || x >= 0)
		return { NAN, NAN };

	// Initialize variables
	mpfr_t xMpfr, low, high;
	mpfr_init2(xMpfr, 53);
	mpfr_init2(low, 53);
	mpfr_init2(high, 53);

	// Convert x to mpfr
	mpfr_set_d(xMpfr, x, MPFR_RNDN);

	// === Compute Bracket ===
	mpfr_t uDown;
	mpfr_init2(uDown, 53);
	mpfr_neg(uDown, xMpfr, MPFR_RNDN);
	mpfr_log(uDown, uDown, MPFR_RNDU);
	mpfr_neg(uDown, uDown, MPFR_RNDN);
	mpfr_sub_ui(uDown, uDown, 1, MPFR_RNDD);

	mpfr_t uUp;
	mpfr_init2(uUp, 53);
	mpfr_neg(uUp, xMpfr, MPFR_RNDN);
	mpfr_log(uUp, uUp, MPFR_RNDD);
	mpfr_neg(uUp, uUp, MPFR_RNDN);
	mpfr_sub_ui(uUp, uUp, 1, MPFR_RNDU);

	// low = -1 - (sqrt(u * 2) + u);
	mpfr_mul_2ui(low, uUp, 1, MPFR_RNDU);
	mpfr_sqrt(low, low, MPFR_RNDU);
	mpfr_add(low, low, uUp, MPFR_RNDU);
	mpfr_si_sub(low, -1, low, MPFR_RNDD);

	// high = -1 - (sqrt(u * 2) + u * 2 / 3)
	mpfr_mul_2ui(high, uDown, 1, MPFR_RNDD);
	mpfr_sqrt(high, high, MPFR_RNDD);
	mpfr_mul_2ui(uDown, uDown, 1, MPFR_RNDD);
	mpfr_div_ui(uDown, uDown, 3, MPFR_RNDD);
	mpfr_add(high, high, uDown, MPFR_RNDD);
	mpfr_si_sub(high, -1, high, MPFR_RNDU);

	mpfr_clear(uDown);
	mpfr_clear(uUp);

	// === Bisection ===
	auto ret = Bisection(x, low, high, false, Wm1Precision);
	assert(ret.sup == std::nextafter(ret.inf, INFINITY));

	mpfr_clear(xMpfr);
	mpfr_clear(low);
	mpfr_clear(high);

	return ret;
}