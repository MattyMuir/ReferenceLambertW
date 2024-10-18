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
	std::uniform_real_distribution<double> dist{ -0.3, 0 };
	for (size_t i = 0; i < arrSize; i++)
		data.push_back(dist(gen));

	// With repeated inits
	TIMER(withInits);
	double _0 = 0;
	for (double d : data)
		_0 += ReferenceWm1(d).inf;
	STOP_LOG(withInits);

	// Without repeated inits
	ReferenceW evaluator;

	TIMER(withoutInits);
	double _1 = 0;
	for (double d : data)
		_1 += evaluator.Wm1(d).inf;
	STOP_LOG(withoutInits);

	// Make sure evaluations are not optimized away (and the same)
	std::cout << std::setprecision(20);
	std::cout << _0 << '\n' << _1 << '\n';
}

int main()
{
	SpeedTest(10000);
}