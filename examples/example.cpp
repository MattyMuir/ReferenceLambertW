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
	std::uniform_real_distribution<double> dist{ 10, 100 };
	for (size_t i = 0; i < arrSize; i++)
		data.push_back(dist(gen));

	// Without repeated inits
	ReferenceW evaluator;

	TIMER(currentVersion);
	double _ = 0;
	for (double d : data)
		_ += evaluator.W0(d).inf;
	STOP_LOG(currentVersion);

	// Make sure evaluations are not optimized away
	std::cout << std::setprecision(20);
	std::cout << _ << '\n';
}

ReferenceW::Sign GetPointSign(double x, double w)
{
	static ReferenceW evaluator;

	// Get sign
	ReferenceW::Sign sign = evaluator.GetMidpointSign(x, w, false);
	if (sign == ReferenceW::Sign::Inconclusive)
		sign = evaluator.GetMidpointSign(x, w, true);

	if (sign == ReferenceW::Sign::Inconclusive)
		throw;

	return sign;
}

void Test()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ 10, 100 };

	ReferenceW evaluator;
	for (;;)
	{
		double x = dist(gen);
		auto [inf, sup] = evaluator.W0(x);

		if (sup != inf && sup != std::nextafter(inf, INFINITY))
			std::cerr << std::format("Bracket too wide! x: {}\n", x);

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
	double x = -1.2885767471783089e-274;

	ReferenceW evaluator;
	evaluator.W0(x);
}