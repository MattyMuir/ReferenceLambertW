#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

#define TIMER_NPRINT
#include "Timer.h"

#include "ReferenceLambertW.h"
#include "ReferenceW.h"

void SpeedTest(size_t arrSize)
{
	// Prepare test data
	std::vector<double> data;
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ 10, 100 };
	for (size_t i = 0; i < arrSize; i++)
		data.push_back(dist(gen));

	// With repeated inits
	TIMER(firstVersion);
	double _0 = 0;
	for (double d : data)
		_0 += ReferenceW0(d).inf;
	STOP_LOG(firstVersion);

	// Without repeated inits
	ReferenceW evaluator;

	TIMER(currentVersion);
	double _1 = 0;
	for (double d : data)
		_1 += evaluator.W0(d).inf;
	STOP_LOG(currentVersion);

	// Make sure evaluations are not optimized away (and the same)
	std::cout << std::setprecision(20);
	std::cout << _0 << '\n' << _1 << '\n';
}

void EqualityTest()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ -0.3678794411714423, 1 };

	ReferenceW evaluator;
	for (;;)
	{
		double x = dist(gen);

		Interval slow = ReferenceW0(x);
		Interval fast = evaluator.W0(x);

		if (slow.inf != fast.inf)
		{
			std::cout << std::setprecision(20);
			std::cout << x << '\n';
		}
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
	SpeedTest(10000);
}