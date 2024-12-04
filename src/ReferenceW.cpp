#include "../include/config.h"
#include "ReferenceW.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cfenv>

#include <iostream>
#include <iomanip>
#include <numeric>

#define SLEEF_STATIC_LIBS
#include <sleef.h>

#include <flint/arb.h>

#include "rndutil.h"
#include "halley.h"

static constexpr double EM_UP = -0.3678794411714423; // (-1/e) rounded towards +Inf
static constexpr double E2_DOWN = 5.43656365691809; // e*2 rounded towards -Inf
static constexpr double E2_UP = 5.436563656918091; // e*2 rounded towards +Inf

ReferenceW::ReferenceW()
{
	arb_init(xArb);
	arb_init(mArb);
	arb_init(yArb);
}

ReferenceW::~ReferenceW()
{
	arb_clear(xArb);
	arb_clear(mArb);
	arb_clear(yArb);
}

Interval ReferenceW::W0(double x)
{
#if REFERENCEW_STATS
	numEvals++;
#endif

	// Edge cases
	if (x < EM_UP)
		return { NAN, NAN };
	if (x == INFINITY)
		return { DBL_MAX, INFINITY };

	// Save current rounding mode
	int initialRnd = fegetround();

	// === Compute Bracket ===
	double high;
	if (x > 3.0)
	{
		// high = ln(x)
		high = Sleef_log_u10(x);
		high = std::nextafter(high, INFINITY);
	}
	else if (x > -0.2875)
	{
		// high = ln(1 + x)
		high = Sleef_log1p_u10(x);
		high = std::nextafter(high, INFINITY);
	}
	else
	{
		// high = -1 + sqrt(2ex + 2)
		high = mul(x, E2_DOWN, FE_UPWARD);
		high = add(high, 2, FE_UPWARD);
		high = sqrt(high, FE_UPWARD);
		high = sub(high, 1, FE_UPWARD);
	}

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
	else if (x > -0.3019)
	{
		// low = x * (1 - x * 2.4)
		double xt = mul(x, 2.4, FE_DOWNWARD);
		low = sub(1, xt, FE_UPWARD);
		low = mul(low, x, FE_DOWNWARD);
	}
	else
	{
		// low = -1 + sqrt(2ex + 2) - 1/3 * (2ex + 2)
		double reta = mul(x, E2_UP, FE_DOWNWARD);
		reta = add(reta, 2, FE_DOWNWARD);
		reta = sqrt(reta, FE_DOWNWARD);

		double eta = mul(x, E2_DOWN, FE_UPWARD);
		eta = add(eta, 2, FE_UPWARD);
		eta = div(eta, 3, FE_UPWARD);

		low = sub(reta, eta, FE_DOWNWARD);
		low = add(low, -1, FE_DOWNWARD);

		// Clamp low above -1
		if (low < -1) low = -1;
	}

	// === Halley Iterations ===
	while (low != 0)
	{
		double newLow = HalleyW0(x, low, false);
		bool stop = (abs((newLow - low) / low) < 1e-7);
		low = newLow;
		if (stop) break;
	}

	while (high != 0)
	{
		double newHigh = HalleyW0(x, high, true);
		bool stop = (abs((newHigh - high) / high) < 1e-7);
		high = newHigh;
		if (stop) break;
	}

	// === Bisection ===
	auto ret = Bisection(x, low, high, true);
	assert(ret.inf == ret.sup || ret.sup == std::nextafter(ret.inf, INFINITY));

	// Restore rounding mode
	fesetround(initialRnd);

	return ret;
}

Interval ReferenceW::Wm1(double x)
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
	auto [low, high] = Wm1Bracket(x);

	// === Bisection ===
	auto ret = Bisection(x, low, high, false);
	assert(ret.inf == ret.sup || ret.sup == std::nextafter(ret.inf, INFINITY));

	// Restore rounding mode
	fesetround(initialRnd);

	return ret;
}

#if REFERENCEW_STATS
double ReferenceW::GetHighPrecRate() const
{
	return (double)numHighPrec / totalBisections;
}

size_t ReferenceW::GetMaxBisections() const
{
	return maxBisections;
}

double ReferenceW::GetAvgBisections() const
{
	return (double)totalBisections / numEvals;
}
#endif

static inline double AddEm(double x)
{
	static constexpr double emHigh = 0.36787944117144232160;
	static constexpr double emLow = -1.2428753672788363168e-17;

	return (x + emHigh) + emLow;
}

static inline double NearBranchWm1(double x)
{
	// === Constants ===
	static constexpr double s2e = 2.331643981597124;
	static constexpr double P[] = {
		-0.9999999999999999,
		-1.0000000000001505,
		-0.3333333333112154,
		-0.15277777908701176,
		-0.07962958804769303,
		-0.0445031235579835,
		-0.02597432503406348,
		-0.015728030108091574,
		-0.009031309914783386,
		-0.008702394675700187,
		0.005169843845676331,
		-0.02414898256188974,
		0.03559281100127844,
		-0.044160933247669634,
		0.030269166389388674,
		-0.011279992858844562
	};
	// =================

	double p = sqrt(AddEm(x)) * s2e;
	double w = P[15];
	for (size_t i = 0; i < 15; i++)
		w = w * p + P[14 - i];

	return w;
}

std::pair<double, double> ReferenceW::Wm1Bracket(double x)
{
	// === Constants ===
	static constexpr double CM13_DOWN = -0.33333333333333337;
	static constexpr double C23_DOWN = 0.6666666666666666;
	static constexpr double C23_UP = 0.6666666666666667;
	static constexpr double N = 50;
	static constexpr double EN_DOWN = 5.184705528587072e+21;
	static constexpr double EN_UP = 5.184705528587073e+21;
	static constexpr double P[] = {
		0,
		-5.415413805902706,
		-2.787876451002007,
		-0.4992978139443087
	};
	static constexpr double Q = 5.410664283026123;
	// =================

	double w;
	if (x > -0.318092372804)
	{
		// Initial approximation
		double t = sqrt(-2 - 2 * log(-x));
		w = P[3];
		for (size_t i = 0; i < 3; i++)
			w = w * t + P[2 - i];
		w = w / (t + Q) - 1.0;

		// Fritsch Iteration
		double zn;
		if (x > -1e-300)
			zn = log((x * 4611686018427387904.0) / w) - 42.975125194716609184 - w;
		else
			zn = log(x / w) - w;
		double temp = 1.0 + w;
		double temp2 = temp + (2.0 / 3.0) * zn;
		temp2 = 2.0 * temp * temp2;
		w = w * (1.0 + (zn / temp) * (temp2 - zn) / (temp2 - 2.0 * zn));
	}
	else
		w = NearBranchWm1(x);

	// Derivative Bound
	fesetround(FE_TONEAREST);
	double logUp = std::nextafter(Sleef_log_u10(-x), INFINITY);
	double rtDown = sqrt(sub(-2, mul(logUp, 2, FE_UPWARD), FE_DOWNWARD), FE_DOWNWARD);
	double denom = add(sub(C23_UP, rtDown, FE_UPWARD), mul(logUp, C23_DOWN, FE_UPWARD), FE_UPWARD);
	double d = sub(1, div(1.0, denom, FE_DOWNWARD), FE_UPWARD);

	// Del Bound
	double del;
	if (x > -0.00000137095397731)
	{
		auto [expDown, expUp] = ExpUpDown(w + N);
		double delDown = mul(div(w, mul(x, EN_UP, FE_DOWNWARD), FE_DOWNWARD), expDown, FE_DOWNWARD);
		double delUp = mul(div(w, mul(x, EN_DOWN, FE_UPWARD), FE_UPWARD), expUp, FE_UPWARD);
		del = std::max(abs(delDown - 1), abs(delUp - 1));
	}
	else
	{
		arb_set_d(mArb, w);
		arb_set_d(xArb, x);
		arb_exp(yArb, mArb, 100);
		arb_mul(yArb, yArb, mArb, 100);
		arb_sub(yArb, yArb, xArb, 100);
		arb_div(yArb, yArb, xArb, 100);
		arb_abs(yArb, yArb);
		arf_t delArf;
		arf_init(delArf);
		arb_get_ubound_arf(delArf, yArb, 100);
		del = arf_get_d(delArf, ARF_RND_UP);
		arf_clear(delArf);
	}

	// Compute final error
	double err = mul(d, del, FE_UPWARD);
	double low = sub(w, err, FE_DOWNWARD);
	double high = add(w, err, FE_UPWARD);
	high = std::min(high, -1.0);

	return { low, high };
}

Sign ReferenceW::GetMidpointSign(double x, double midpoint, bool useHighPrec)
{
#if REFERENCEW_STATS
	if (useHighPrec)
		numHighPrec++;
#endif

#if 0
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

	return Sign::Inconclusive;
#else
	slong prec = useHighPrec ? 150 : 90;

	arb_set_d(xArb, x);
	arb_set_d(mArb, midpoint);
	arb_exp(yArb, mArb, prec);
	arb_mul(yArb, yArb, mArb, prec);
	arb_sub(yArb, yArb, xArb, prec);

	bool isPos = arb_is_nonnegative(yArb);
	bool isNeg = arb_is_nonpositive(yArb);

	if (isPos)
		return Sign::Positive;
	if (isNeg)
		return Sign::Negative;

	return Sign::Inconclusive;
#endif
}

Interval ReferenceW::Bisection(double x, double low, double high, bool increasing)
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

#if REFERENCEW_STATS
	maxBisections = std::max(maxBisections, b);
	totalBisections += b;
#endif

	return { low, high };
}