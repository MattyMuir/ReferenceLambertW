#include "../include/config.h"
#include "ReferenceWf.h"

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
	arb_init(xArb);
	arb_init(mArb);
	arb_init(yArb);
}

ReferenceWf::~ReferenceWf()
{
	arb_clear(xArb);
	arb_clear(mArb);
	arb_clear(yArb);
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
	if (x == 0)
		return { 0, 0 };

	// Save current rounding mode
	int initialRnd = fegetround();

	// === Compute Bracket ===
	auto [low, high] = W0Bracket(x);

	// === Bisection ===
	auto ret = Bisection(x, low, high, true);
	if (!(ret.inf == ret.sup || ret.sup == std::nextafter(ret.inf, INFINITY)))
	{
		std::cerr << std::format("Bracket too wide x: {}\n", x);
		std::terminate();
	}

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
	auto [low, high] = Wm1Bracket(x);

	// === Bisection ===
	auto ret = Bisection(x, low, high, false);
	if (!(ret.inf == ret.sup || ret.sup == std::nextafter(ret.inf, INFINITY)))
	{
		std::cerr << std::format("Bracket too wide x: {}\n", x);
		std::terminate();
	}

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

static inline float AddEm(float x)
{
	static constexpr float emHigh = 0.36787945f;
	static constexpr float emLow = -9.149756e-09f;
	return (x + emHigh) + emLow;
}

static inline float FirstApproxW0(float x)
{
	static constexpr double P[] = {
		0,
		165.51561672164559,
		1104.9153130867758,
		2632.284078577963,
		2689.464120405435,
		1121.2923665114324,
		153.3374641092571,
		4.077322829553558
	};

	static constexpr double Q[] = {
		165.51561558818844,
		1270.4310030077481,
		3654.442208397931,
		4879.631928655197,
		3045.0058891120098,
		794.8712729472717,
		67.22857835896016,
		1
	};

	double numer = P[7];
	for (size_t i = 0; i < 7; i++)
		numer = numer * x + P[6 - i];

	double denom = Q[7];
	for (size_t i = 0; i < 7; i++)
		denom = denom * x + Q[6 - i];

	return numer / denom;
}

static inline float SecondApproxW0(float x)
{
	static constexpr double P[] = {
		245182.20097823755,
		280243.5212428723,
		142843.813324628,
		40353.72076097795,
		5776.914448840662,
		184.83613670644033,
		0.9984483567344636
	};

	static constexpr double Q[] = {
		432788.26007218857,
		216948.13159273885,
		58081.26591912717,
		6594.751582203545,
		191.21022696372594,
		1
	};

	double t = log((double)x);

	double numer = P[6];
	for (size_t i = 0; i < 6; i++)
		numer = numer * t + P[5 - i];

	double denom = Q[5];
	for (size_t i = 0; i < 5; i++)
		denom = denom * t + Q[4 - i];

	return numer / denom;
}

static inline float NearBranchW0(float x)
{
	static constexpr double e2 = 5.43656365691809;

	static constexpr double P[] = {
		-0.9999999781289544,
		0.9999966080647236,
		-0.33324531164727067,
		0.15189891604646868,
		-0.07530393941472714,
		0.03290035332102544,
		-0.008369773627101843
	};

	double p = sqrt(e2 * x + 2.0);

	double res = P[6];
	for (size_t i = 0; i < 6; i++)
		res = res * p + P[5 - i];

	return res;
}

std::pair<float, float> ReferenceWf::W0Bracket(float x)
{
	float w = (x < -0.3f) ? NearBranchW0(x) : ((x < 7.38905609893f) ? FirstApproxW0(x) : SecondApproxW0(x));

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
	else
	{
		if (x < -0.01)
		{
			static constexpr double a = -0.1321205588285577; // (2 - e) / 2e rounded towards -Inf
			static constexpr double b = 0.8939534673502061; // sqrt(2)(e - 1) / e rounded towards -Inf
			static constexpr double E2_DOWN = 5.43656365691809;
			static constexpr double E2_UP = 5.436563656918091;

			double etaUp = fma(E2_DOWN, (double)x, 2, FE_UPWARD);
			double etaDown = fma(E2_UP, (double)x, 2, FE_DOWNWARD);
			etaDown = mul(b, sqrt(etaDown, FE_DOWNWARD), FE_DOWNWARD);
			double denom = fma(a, etaUp, etaDown, FE_DOWNWARD);
			d = sub(div(1, denom, FE_UPWARD), 1, FE_UPWARD);
		}
		else
			d = sub(mul(mul((double)x, (double)x, FE_UPWARD), 3, FE_UPWARD), (double)x, FE_UPWARD);
	}

	// Del bound
	auto [expDown, expUp] = ExpUpDown((double)w);
	double delDown = mul(div((double)w, (double)x, FE_DOWNWARD), expDown, FE_DOWNWARD);
	double delUp = mul(div((double)w, (double)x, FE_UPWARD), expUp, FE_UPWARD);
	double del = std::max(abs(delDown - 1), abs(delUp - 1));

	// Compute final error
	float err = mul(d, del, FE_UPWARD);
	float low = sub(w, err, FE_DOWNWARD);
	float high = add(w, err, FE_UPWARD);
	high = std::max(high, -1.0f);

	if (low == 0) low = 0;

	return { low, high };
}

static inline float NearBranchWm1(float x)
{
	static constexpr float s2e = 2.331644f;

	float p = s2e * sqrt(AddEm(x));

	static constexpr float P[] = {
		-1.0000000001291165,
		-0.9999992250595189,
		-0.3340219624089988
	};

	float res = P[2];
	for (size_t i = 0; i < 2; i++)
		res = res * p + P[1 - i];

	return res;
}

static inline float GeneralWm1(float x)
{
	static constexpr double P[] = {
		-2101.555169658076,
		-3413.0457024602106,
		-2345.4071921263444,
		-864.1804177336671,
		-175.99964384176346,
		-17.64071303855079,
		-0.4998769261313046
	};

	static constexpr double Q[] = {
		2101.5551872949245,
		1311.4898275251383,
		333.4030604186147,
		35.228646667156625,
		1
	};

	double t = sqrt(-2 - 2 * log((double)-x));

	double numer = P[6];
	for (size_t i = 0; i < 6; i++)
		numer = numer * t + P[5 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

std::pair<float, float> ReferenceWf::Wm1Bracket(float x)
{
	// === Constants ===
	static constexpr double C23_DOWN = 0.6666666666666666;
	static constexpr double C23_UP = 0.6666666666666667;
	// =================

	float w = (x < -0.367877785718f) ? NearBranchWm1(x) : GeneralWm1(x);

	// Derivative Bound
	double logUp = std::nextafter(Sleef_log_u10(-x), INFINITY);
	double rtDown = sqrt(sub(-2, mul(logUp, 2, FE_UPWARD), FE_DOWNWARD), FE_DOWNWARD);
	double denom = add(sub(C23_UP, rtDown, FE_UPWARD), mul(logUp, C23_DOWN, FE_UPWARD), FE_UPWARD);
	double d = sub(1, div(1.0, denom, FE_DOWNWARD), FE_UPWARD);

	// Del bound
	auto [expDown, expUp] = ExpUpDown((double)w);
	double delDown = mul(div((double)w, (double)x, FE_DOWNWARD), expDown, FE_DOWNWARD);
	double delUp = mul(div((double)w, (double)x, FE_UPWARD), expUp, FE_UPWARD);
	double del = std::max(abs(delDown - 1), abs(delUp - 1));

	// Compute final error
	float err = mul(d, del, FE_UPWARD);
	float low = sub(w, err, FE_DOWNWARD);
	float high = add(w, err, FE_UPWARD);
	high = std::min(high, -1.0f);

	return { low, high };
}

Sign ReferenceWf::GetMidpointSign(float x, float midpoint, bool useHighPrec)
{
#if REFERENCEW_STATS
	if (useHighPrec)
		numHighPrec++;
#endif

	if (!useHighPrec)
	{
		if (midpoint >= x)
			return Sign::Positive;

		double m = midpoint;

		// Compute exp
		auto [yLow, yHigh] = ExpUpDown(m);
		if (midpoint < 0)
			std::swap(yLow, yHigh);

		// Compute yLow
		yLow = mul(yLow, m, FE_DOWNWARD);
		yLow = sub(yLow, (double)x, FE_DOWNWARD);

		// Compute yHigh
		yHigh = mul(yHigh, m, FE_UPWARD);
		yHigh = sub(yHigh, (double)x, FE_UPWARD);

		if (yLow >= 0 && yHigh >= 0)
			return Sign::Positive;
		if (yLow <= 0 && yHigh <= 0)
			return Sign::Negative;

		return Sign::Inconclusive;
	}

	// High precision implementation
	arb_set_d(xArb, x);
	arb_set_d(mArb, midpoint);
	arb_exp(yArb, mArb, 70);
	arb_mul(yArb, yArb, mArb, 70);
	arb_sub(yArb, yArb, xArb, 70);

	bool isPos = arb_is_nonnegative(yArb);
	bool isNeg = arb_is_nonpositive(yArb);

	if (isPos)
		return Sign::Positive;
	if (isNeg)
		return Sign::Negative;

	return Sign::Inconclusive;
}

Intervalf ReferenceWf::Bisection(float x, float low, float high, bool increasing)
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
		float m = std::midpoint(low, high);

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