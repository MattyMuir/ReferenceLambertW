#include "../include/config.h"
#include "ReferenceW.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cfenv>

#include <iostream>
#include <format>
#include <numeric>

#define SLEEF_STATIC_LIBS
#include <sleef.h>

#include <arb.h>

#include "rndutil.h"
#include "halley.h"

static constexpr double EM_UP = -0.3678794411714423; // (-1/e) rounded towards +Inf

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
	if (x == 0)
		return { 0, 0 };

	// Save current rounding mode
	int initialRnd = fegetround();

	// === Compute Bracket ===
	auto [low, high] = W0Bracket(x);

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

static inline double NearBranchW0(double x)
{
	static constexpr double s2e = 2.331643981597124;
	static constexpr double P[] = {
		-1.00000000000000000000,
		0.99999999999998689937,
		-0.33333333333171155655,
		0.15277777769847986078,
		-0.07962962759798784818,
		0.04450228328389740917,
		-0.02598439214142129680,
		0.01563333375832150554,
		-0.00960508856297833703,
		0.00596982547465134492,
		-0.00368441824865070513,
		0.00216878673408957843,
		-0.00113330227139719539,
		0.00047252681627728467,
		-0.00013420111092875102,
		0.00001887878365359131,
	};

	double p = sqrt(AddEm(x)) * s2e;

	double value = P[15];
	for (size_t i = 0; i < 15; i++)
		value = value * p + P[14 - i];

	return value;
}

static inline double FirstW0Approx(double x)
{
	if (abs(x) < 1e-4)
		return x;

	static constexpr double P[] = {
		0,
		30.580056454638136,
		83.95836185597197,
		46.16620637664877,
		3.4636816277252214
	};

	static constexpr double Q[] = {
		30.578403642151667,
		114.49011569793561,
		114.80618615998705,
		28.635096582884064,
		1
	};

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * x + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * x + Q[3 - i];

	return numer / denom;
}

double SecondW0Approx(double x)
{
	static constexpr double P[] = {
		64312.7454007891,
		43264.12227598657,
		20243.65384336377,
		453.17656235798086,
		1.0000432316050645
	};
	static constexpr double Q[] = {
		104342.57917932322,
		22499.368605590193,
		460.93750724715477,
		1
	};

	double lx = log(x);

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * lx + P[3 - i];

	double denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * lx + Q[2 - i];

	return numer / denom;
}

std::pair<double, double> ReferenceW::W0Bracket(double x)
{
	// Initial approximation
	double w;
	fesetround(FE_TONEAREST);
	if (x < -0.28)
		w = NearBranchW0(x);
	else
	{
		w = (x < 7.34) ? FirstW0Approx(x) : SecondW0Approx(x);

		// Fritsch Iteration
		double zn = log(x / w) - w;
		double temp = 1.0 + w;
		double temp2 = temp + (2.0 / 3.0) * zn;
		temp2 = 2.0 * temp * temp2;
		w = w * (1.0 + (zn / temp) * (temp2 - zn) / (temp2 - 2.0 * zn));
	}

	// Derivative Bound
	double d = x;
	if (x > 0)
	{
		if (x > 0.01)
		{
			double logUp = Sleef_log1p_u10(x);
			logUp = std::nextafter(logUp, INFINITY);
			d = sub(1, div(1, add(1, logUp, FE_UPWARD), FE_DOWNWARD), FE_UPWARD);
		}
	}
	else if (x < 0)
	{
		if (x < -0.01)
		{
			static constexpr double a = -0.1321205588285577; // (2 - e) / 2e rounded towards -Inf
			static constexpr double b = 0.8939534673502061; // sqrt(2)(e - 1) / e rounded towards -Inf
			static constexpr double E2_DOWN = 5.43656365691809;
			static constexpr double E2_UP = 5.436563656918091;

			double etaUp = fma(E2_DOWN, x, 2, FE_UPWARD);
			double etaDown = fma(E2_UP, x, 2, FE_DOWNWARD);
			etaDown = mul(b, sqrt(etaDown, FE_DOWNWARD), FE_DOWNWARD);
			double denom = fma(a, etaUp, etaDown, FE_DOWNWARD);
			d = sub(div(1, denom, FE_UPWARD), 1, FE_UPWARD);
		}
		else
			d = sub(mul(mul(x, x, FE_UPWARD), 3, FE_UPWARD), x, FE_UPWARD);
	}

	// Del Bound
	double del;
	if (x > 4.11380962917)
	{
		static constexpr double N = 50;
		static constexpr double EN_DOWN = 5.184705528587072e+21;
		static constexpr double EN_UP = 5.184705528587073e+21;

		if (x > 1)
		{
			auto [expDown, expUp] = ExpUpDown(w);
			double delDown = mul(div(w, x, FE_DOWNWARD), expDown, FE_DOWNWARD);
			double delUp = mul(div(w, x, FE_UPWARD), expUp, FE_UPWARD);
			del = std::max(abs(delDown - 1), abs(delUp - 1));
		}
		else
		{
			double expDown = add(w, N, FE_DOWNWARD);
			double expUp = add(w, N, FE_UPWARD);
			fesetround(FE_TONEAREST);
			expDown = std::nextafter(Sleef_exp_u10(expDown), -INFINITY);
			expUp = std::nextafter(Sleef_exp_u10(expUp), INFINITY);
			double delDown = mul(div(w, mul(x, EN_UP, FE_DOWNWARD), FE_DOWNWARD), expDown, FE_DOWNWARD);
			double delUp = mul(div(w, mul(x, EN_DOWN, FE_UPWARD), FE_UPWARD), expUp, FE_UPWARD);
			del = std::max(abs(delDown - 1), abs(delUp - 1));
		}
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
	high = std::max(high, -1.0);

	return { low, high };
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
	fesetround(FE_TONEAREST);
	if (x < -0.318092372804)
		w = NearBranchWm1(x);
	else
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


	// Derivative Bound
	double logUp = std::nextafter(Sleef_log_u10(-x), INFINITY);
	double rtDown = sqrt(sub(-2, mul(logUp, 2, FE_UPWARD), FE_DOWNWARD), FE_DOWNWARD);
	double denom = add(sub(C23_UP, rtDown, FE_UPWARD), mul(logUp, C23_DOWN, FE_UPWARD), FE_UPWARD);
	double d = sub(1, div(1.0, denom, FE_DOWNWARD), FE_UPWARD);

	// Del Bound
	double del;
	if (x > -0.00000137095397731)
	{
		double expDown = add(w, N, FE_DOWNWARD);
		double expUp = add(w, N, FE_UPWARD);
		fesetround(FE_TONEAREST);
		expDown = std::nextafter(Sleef_exp_u10(expDown), -INFINITY);
		expUp = std::nextafter(Sleef_exp_u10(expUp), INFINITY);
		double delDown = mul(div(w, mul(x, EN_UP, FE_DOWNWARD), FE_DOWNWARD), expDown, FE_DOWNWARD);
		delDown = sub(delDown, 1, FE_DOWNWARD);
		double delUp = mul(div(w, mul(x, EN_DOWN, FE_UPWARD), FE_UPWARD), expUp, FE_UPWARD);
		delUp = sub(delUp, 1, FE_UPWARD);
		del = std::max(abs(delDown), abs(delUp));
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

	if (midpoint >= x)
		return Sign::Positive;
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
}

Interval ReferenceW::Bisection(double x, double low, double high, bool increasing)
{
#if REFERENCEW_STATS
	size_t b = 0;
#endif

	fesetround(FE_TONEAREST);
	for (;;)
	{
#if REFERENCEW_STATS
		b++;
#endif

		if (high <= std::nextafter(low, INFINITY))
			break; // Bracket cannot be narrowed any further

		// m = (low + high) / 2
		double m = std::midpoint(low, high);

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