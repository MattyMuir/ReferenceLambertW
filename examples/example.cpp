#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

#define TIMER_NPRINT
#include "Timer.h"

#include <ReferenceLambertW.h>

template <typename Ty>
void SpeedTest(size_t arrSize)
{
	// Prepare test data
	std::vector<Ty> data;
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<Ty> dist{ -0.35, -0.2 };
	for (size_t i = 0; i < arrSize; i++)
		data.push_back(dist(gen));

	// Create evaluator object
	std::conditional_t<std::is_same_v<Ty, float>, ReferenceWf, ReferenceW> evaluator;

	TIMER(t);
	Ty _ = 0;
	for (Ty d : data)
		_ += evaluator.Wm1(d).inf;
	STOP_LOG(t);

	// Make sure evaluations are not optimized away
	std::cout << std::setprecision(20);
	std::cout << _ << '\n';
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
	double x = -0.36;

	ReferenceW evaluator;
	auto[inf, sup] = evaluator.Wm1(x);
	std::cout << std::format("{} - {}", inf, sup);
}