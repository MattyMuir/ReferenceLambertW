#include <iostream>
#include <string>
#include <charconv>
#include <random>
#include <format>

#include <mpfr.h>
#include <ReferenceLambertW.h>

#include "ReciprocalDistributionEx.h"

#define ERROR(msg) { std::cout << msg << '\n'; return 1; }

template <typename Ty>
Ty GetEmUp()
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
int RunTest(int64_t branch)
{
	Ty low = GetEmUp<Ty>();
	Ty high = (branch == 0) ? INFINITY : 0;

	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx<Ty> dist{ low, high, false };

	// Construct evaluators
	ReferenceW evaluatord;
	ReferenceWf evaluatorf;

	for (size_t i = 0; i < 500'000; i++)
	{
		Ty x = dist(gen);
		std::cout << x << '\n';
		
		if constexpr (std::is_same_v<Ty, float>)
		{
			Intervalf res;
			if (branch == 0)
				res = evaluatorf.W0(x);
			else
				res = evaluatorf.Wm1(x);

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
		}
		else
		{
			Interval res;
			if (branch == 0)
				res = evaluatord.W0(x);
			else
				res = evaluatord.Wm1(x);

			if (res.inf != res.sup && res.sup != std::nextafter(res.inf, INFINITY))
			{
				std::cerr << std::format("Too wide x: {}\n", x);
				return 1;
			}

			if (WexpwIsPositive(res.inf, x) == WexpwIsPositive(res.sup, x))
			{
				std::cerr << std::format("x: {}\n", x);
				return 1;
			}
		}
	}

	return 0;
}

int main(int argc, char** argv)
{
#if 0
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
#else
	size_t testIdx = 0;
#endif

	switch (testIdx)
	{
	case 0: return RunTest<float>(0);
	case 1: return RunTest<float>(-1);
	case 2: return RunTest<double>(0);
	case 3: return RunTest<double>(-1);
	default: ERROR("Invalid test index");
	}
}