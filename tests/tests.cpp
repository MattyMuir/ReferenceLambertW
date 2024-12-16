#include <iostream>
#include <string>
#include <charconv>
#include <random>
#include <format>
#include <functional>

#include <mpfr.h>
#include <ReferenceLambertW.h>

#include "ReciprocalDistributionEx.h"

#define ERROR(msg) { std::cerr << msg << '\n'; return 1; }

template <typename Ty>
consteval Ty GetEmUp()
{
	if constexpr (std::is_same_v<Ty, float>)
		return -0.36787942f;
	else
		return -0.3678794411714423;
}

template <typename Ty>
mpfr_prec_t GetPrec()
{
	if constexpr (std::is_same_v<Ty, float>)
		return 24;
	else
		return 53;
}

template <typename Ty>
bool WexpwIsPositive(Ty w, Ty x)
{
	mpfr_t wMpfr, yLow, yHigh;
	mpfr_init2(wMpfr, GetPrec<Ty>());
	mpfr_init2(yLow, 150);
	mpfr_init2(yHigh, 150);

	// Convert w to mpfr
	if constexpr (std::is_same_v<Ty, float>)
		mpfr_set_flt(wMpfr, w, MPFR_RNDN);
	else
		mpfr_set_d(wMpfr, w, MPFR_RNDN);

	// Compute exp
	mpfr_exp(yLow, wMpfr, MPFR_RNDD);
	mpfr_set(yHigh, yLow, MPFR_RNDN);
	mpfr_nextabove(yHigh);
	if (mpfr_cmp_ui(wMpfr, 0) < 0)
		mpfr_swap(yLow, yHigh);

	// Compute yLow
	mpfr_mul(yLow, yLow, wMpfr, MPFR_RNDD);
	mpfr_sub_d(yLow, yLow, x, MPFR_RNDD);

	// Compute yHigh
	mpfr_mul(yHigh, yHigh, wMpfr, MPFR_RNDU);
	mpfr_sub_d(yHigh, yHigh, x, MPFR_RNDU);

	int lowCmp = mpfr_cmp_ui(yLow, 0);
	int highCmp = mpfr_cmp_ui(yHigh, 0);

	mpfr_clear(wMpfr);
	mpfr_clear(yLow);
	mpfr_clear(yHigh);

	if (lowCmp >= 0 && highCmp >= 0)
		return true;
	if (lowCmp <= 0 && highCmp <= 0)
		return false;

	std::cerr << "Inconclusive sign!\n";
	throw;
}

template <typename Ty>
int TestPoint(Ty x, const std::conditional_t<std::is_same_v<Ty, float>, Intervalf, Interval>& res)
{
	if (res.inf != res.sup && res.sup != std::nextafter(res.inf, INFINITY))
	{
		std::cerr << std::format("Too wide x: {}\n", x);
		return 1;
	}

	if (WexpwIsPositive(res.inf, x) == WexpwIsPositive(res.sup, x))
	{
		std::cerr << std::format("Incorrect x: {}\n", x);
		return 1;
	}

	return 0;
}

template <typename Ty>
int RunTest(int64_t branch, const std::function<Ty()>& rand)
{
	// === Parameters ===
	static constexpr size_t Num = 500'000;
	// ==================

	// Construct evaluators
	std::conditional_t<std::is_same_v<Ty, float>, ReferenceWf, ReferenceW> evaluator;

	for (size_t i = 0; i < Num; i++)
	{
		Ty x = rand();
		
		if constexpr (std::is_same_v<Ty, float>)
		{
			Intervalf res;
			if (branch == 0)
				res = evaluator.W0(x);
			else
				res = evaluator.Wm1(x);

			if (TestPoint(x, res)) return 1;
		}
		else
		{
			Interval res;
			if (branch == 0)
				res = evaluator.W0(x);
			else
				res = evaluator.Wm1(x);

			if (TestPoint(x, res)) return 1;
		}
	}

	return 0;
}

template <typename Ty>
int RunTest(int64_t branch)
{
	static std::mt19937_64 gen{ std::random_device{}() };

	// ReciprocalDist test
	{
		Ty low = GetEmUp<Ty>();
		Ty high = (branch == 0) ? INFINITY : 0;

		ReciprocalDistributionEx<Ty> dist{ low, high, false };
		if (RunTest<Ty>(branch, [&]() { return dist(gen); }))
			return 1;
	}

	// Near branch test
	{
		Ty low = std::is_same_v<Ty, float> ? -15 : -35;
		Ty high = -1;

		std::uniform_real_distribution<Ty> dist{ low, high };
		if (RunTest<Ty>(branch, [&]() { return GetEmUp<Ty>() + exp(dist(gen)); }))
			return 1;
	}

	// Zero test
	if (branch == 0)
	{
		std::conditional_t<std::is_same_v<Ty, float>, ReferenceWf, ReferenceW> evaluator;
		auto [inf, sup] = evaluator.W0(0);
		if (inf != 0 || sup != 0)
		{
			std::cerr << "Failed zero test!\n";
			return 1;
		}
	}

	return 0;
}

template <typename Ty>
int ExhaustiveTest(int64_t branch)
{
	std::conditional_t<std::is_same_v<Ty, float>, ReferenceWf, ReferenceW> evaluator;

	Ty start = GetEmUp<Ty>();
	Ty end = (branch == 0) ? INFINITY : 0;

	Ty x = start;
	for (size_t i = 0; x < end; i++, x = std::nextafter(x, INFINITY))
	{
		decltype(evaluator.W0(Ty{})) res;
		if (branch == 0)
			res = evaluator.W0(x);
		else
			res = evaluator.Wm1(x);

		if (TestPoint(x, res)) return 1;
	}

	return 0;
}

int main(int argc, char** argv)
{
	// Check number of arguments is correct
	if (argc != 2)
		ERROR("Test must have exactly one extra argument");

	// Get second argument
	std::string arg{ argv[1] };

	// Convert to integer
	size_t testIdx;
	auto convRes = std::from_chars(arg.data(), arg.data() + arg.size(), testIdx);
	if (convRes.ec != std::errc())
		ERROR("Test index could not be parsed");

	switch (testIdx)
	{
	case 0: return RunTest<float>(0);
	case 1: return RunTest<float>(-1);
	case 2: return RunTest<double>(0);
	case 3: return RunTest<double>(-1);
	case 4: return ExhaustiveTest<float>(0);
	case 5: return ExhaustiveTest<float>(-1);
	default: ERROR("Invalid test index");
	}
}