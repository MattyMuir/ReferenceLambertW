#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

#define TIMER_NPRINT
#include "Timer.h"

#include "ReferenceW.h"

void SpeedTest(size_t arrSize)
{
	// Prepare test data
	std::vector<double> data;
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ -0.3, 0 };
	for (size_t i = 0; i < arrSize; i++)
		data.push_back(dist(gen));

	// Without repeated inits
	ReferenceW evaluator;

	TIMER(currentVersion);
	double _ = 0;
	for (double d : data)
		_ += evaluator.Wm1(d).inf;
	STOP_LOG(currentVersion);

	// Make sure evaluations are not optimized away
	std::cout << std::setprecision(20);
	std::cout << _ << '\n';
}

ReferenceW::Sign GetPointSign(double x, double w)
{
	static ReferenceW evaluator;

	// Convert w to MPFR
	mpfr_t wMpfr;
	mpfr_init2(wMpfr, 53);
	mpfr_set_d(wMpfr, w, MPFR_RNDN);

	// Get sign
	ReferenceW::Sign sign = evaluator.GetMidpointSign(x, wMpfr, false);
	if (sign == ReferenceW::Sign::Inconclusive)
		sign = evaluator.GetMidpointSign(x, wMpfr, true);

	// Clear MPFR variable
	mpfr_clear(wMpfr);

	if (sign == ReferenceW::Sign::Inconclusive)
		throw;

	return sign;
}

void Test()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ -1e-310, -1e-315 };

	ReferenceW evaluator;
	for (;;)
	{
		double x = dist(gen);
		auto [inf, sup] = evaluator.Wm1(x);

		auto infSign = GetPointSign(x, inf);
		auto supSign = GetPointSign(x, sup);
		if (infSign == supSign)
			std::cerr << std::format("Error! x: {}\n", x);
	}
}

void Profiling(size_t arrSize, size_t numIter)
{
	// Prepare test data
	std::vector<double> data;
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ 10, 100 };
	for (size_t i = 0; i < arrSize; i++)
		data.push_back(dist(gen));

	ReferenceW evaluator;

	double _ = 0.0;
	for (size_t i = 0; i < numIter; i++)
		for (double d : data)
			_ += evaluator.W0(d).inf;

	std::cout << _;
}

int main()
{
	std::vector<double> xs;
	for (double x = 0.24; x < 709; x += 0.5)
		xs.push_back(x);

	std::vector<double> ys = xs;
	ReferenceW evaluator;
	std::transform(ys.begin(), ys.end(), ys.begin(), [&](double v) { return evaluator.Wm1(-exp(-v - 1)).inf; });

	for (double x : xs)
		std::cout << std::format("{}\n", x);

	std::cout << "\n\n\n\n\n";

	for (double y : ys)
		std::cout << std::format("{}\n", y);
}