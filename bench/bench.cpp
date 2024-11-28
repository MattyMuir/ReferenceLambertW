#include <iostream>
#include <vector>
#include <fstream>
#include <random>

#include <ReferenceLambertW.h>

#define TIMER_NPRINT
#include "Timer.h"

// === Bench Config ===
#define BRANCH W0
using BenchTy = double;
// ====================

#define DOMAP(x) ExpMap##x
#define MAP(x) DOMAP(x)

template <typename Ty>
using Function1D = Ty(*)(Ty);

template <typename Ty>
double RunBench(Ty min, Ty max, Function1D<Ty> map, size_t num)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<Ty> dist{ min, max };

	// Prepare data
	std::vector<Ty> data;
	data.reserve(num);
	for (size_t i = 0; i < num; i++)
		data.push_back(map(dist(gen)));

	// Create evaluator
	static std::conditional_t<std::is_same_v<Ty, float>, ReferenceWf, ReferenceW> evaluator;

	// Run timing
	Ty _ = 0;
	Timer t;
	for (Ty d : data)
		_ += evaluator.BRANCH(d).inf;
	t.Stop();

	return t.GetSeconds();
}

template <typename Ty>
std::tuple<double, size_t, double> RunStats(Ty min, Ty max, Function1D<Ty> map, size_t num)
{
#if REFERENCEW_STATS
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<Ty> dist{ min, max };

	// Create evaluator
	std::conditional_t<std::is_same_v<Ty, float>, ReferenceWf, ReferenceW> evaluator;

	// Track stats
	for (size_t i = 0; i < num; i++)
		evaluator.BRANCH(map(dist(gen)));

	return { evaluator.GetHighPrecRate(), evaluator.GetMaxBisections(), evaluator.GetAvgBisections() };
#else
	return { 0.0, 0, 0.0 };
#endif
}

float ExpMapW0(float x)
{
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP + exp(x);
}

float ExpMapWm1(float x)
{
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP / (1 + exp(x));
}

double ExpMapW0(double x)
{
	static constexpr double EM_UP = -0.3678794411714423;
	return EM_UP + exp(x);
}

double ExpMapWm1(double x)
{
	static constexpr double EM_UP = -0.3678794411714423;
	return EM_UP / (1 + exp(x));
}

void Bench()
{
	// === Parameters ===
	static constexpr size_t Num = 10'000;
	static constexpr size_t Repeats = 30;
	BenchTy binMin = -10;
	BenchTy binMax = 400;
	BenchTy binWidth = 1;
	// ==================

	std::ofstream file{ "bench.csv" };

	file << "Min,Max,Time\n";
	for (BenchTy min = binMin; min < binMax; min += binWidth)
	{
		BenchTy max = min + binWidth;
		double time = 0.0;
		for (size_t i = 0; i < Repeats; i++)
			time += RunBench(min, max, MAP(BRANCH), 100);
		time /= Repeats;

		file << std::format("{:.3f},{:.3f},{:.10f}\n", min, max, time);
		std::cout << min << " - " << max << '\n';
	}
}

void Stats()
{
#if REFERENCEW_STATS
	// === Parameters ===
	static constexpr size_t Num = 1'000;
	BenchTy binMin = -10;
	BenchTy binMax = 400;
	BenchTy binWidth = 1;
	// ==================

	std::ofstream file{ "stats.csv" };

	file << "Min,Max,HighPrec Rate,Max Bisections,Average Bisections\n";
	for (BenchTy min = binMin; min < binMax; min += binWidth)
	{
		BenchTy max = min + binWidth;
		auto [highPrecRate, maxBisections, avgBisections] = RunStats(min, max, MAP(BRANCH), Num);

		file << std::format("{:.3f},{:.3f},{:.10f},{},{:.10f}\n", min, max, highPrecRate, maxBisections, avgBisections);
		file << std::flush;
		std::cout << min << " - " << max << '\n';
	}
#endif
}

int main()
{
	Bench();
	Stats();
	//double time = RunBench()
}