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

void Test()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ -0.3678794411714423, 0 };

	ReferenceW evaluator;
	for (;;)
	{
		// TODO
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
	SpeedTest(10'000);
}