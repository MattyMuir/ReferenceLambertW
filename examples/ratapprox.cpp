#include <cassert>

#include <iostream>
#include <vector>
#include <format>
#include <random>
#include <thread>

#include <ReferenceLambertW.h>

// === Parameters ===
constexpr size_t pOrder = 7;
constexpr size_t qOrder = 7;
constexpr size_t numCoeffs = pOrder + qOrder + 1;
constexpr double yOffset = 0;

constexpr double min = -0.3;
constexpr double max = 6.9035267829895019531;
constexpr double step = 0.005;
// ==================

struct DerivState
{
	double x, y, numer, denom, signedError;
};

std::vector<double> GetCoefficients(const std::vector<double>& ps, const std::vector<double>& scales)
{
	std::vector<double> coeffs = ps;
	for (size_t i = 0; i < coeffs.size(); i++)
		coeffs[i] *= scales[i];

	return coeffs;
}

static double IntPow(double x, size_t p)
{
	double ret = 1;
	for (size_t i = 0; i < p; i++)
		ret *= x;
	return ret;
}

static double sign(double x)
{
	return (x > 0) ? 1.0 : -1.0;
}

template <typename Ty>
Ty EvaluateRational(Ty x, const std::vector<Ty>& coeffs, Ty* numer_ = nullptr, Ty* denom_ = nullptr)
{
	// Evaluate rational
	Ty numer = coeffs[pOrder];
	for (size_t i = 0; i < pOrder; i++)
		numer = numer * x + coeffs[pOrder - 1 - i];

	Ty denom = 1;
	for (size_t i = 0; i < qOrder; i++)
		denom = denom * x + coeffs[pOrder + qOrder - i];

	// Calculate error
	Ty approx = numer / denom + (Ty)yOffset;

	if (numer_) *numer_ = numer;
	if (denom_) *denom_ = denom;

	return approx;
}

std::pair<double, std::vector<double>> GetMaxError(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& ps, const std::vector<double>& scales)
{
	// Calculate coefficients
	std::vector<double> coeffs = GetCoefficients(ps, scales);

	// Calculate max error
	double maxError = 0;
	DerivState state;
	for (size_t i = 0; i < xs.size(); i++)
	{
		double x = xs[i];
		double numer, denom;
		double approx = EvaluateRational(x, coeffs, &numer, &denom);
		
		double y = ys[i];
		double signedError = (approx - y) / y;
		double error = abs(signedError);
		if (error > maxError)
		{
			maxError = error;
			state = DerivState{ x, y, numer, denom, signedError };
		}
	}

	// Calculate derivatives
	std::vector<double> derivs;
	for (size_t pi = 0; pi < numCoeffs; pi++)
	{
		bool isNumer = (pi <= pOrder);

		double deriv;
		if (isNumer)
			deriv = IntPow(state.x, pi) / (state.y * state.denom) * sign(state.signedError);
		else
			deriv = -IntPow(state.x, pi - pOrder - 1) * state.numer / (state.y * state.denom * state.denom) * sign(state.signedError);

		deriv = (scales[pi] == 0) ? 0 : deriv * scales[pi];
		derivs.push_back(deriv);
	}

	return { maxError, derivs };
}

double W0(double x)
{
	static ReferenceW evaluator;
	return evaluator.W0(x).inf;
}

double Wm1(double x)
{
	static ReferenceW evaluator;
	return evaluator.Wm1(x).inf;
}

double Func(double x)
{
	//return Wm1(-exp(-0.5 * (x * x + 2)));
	//return W0(exp(x));
	return W0(x);
}

float fFunc(float x)
{
	return Func(x);
}

uint32_t ULPDistance(float a, float b)
{
	if (a > b)
		std::swap(a, b);
	uint32_t aPunn = std::bit_cast<uint32_t>(a);
	uint32_t bPunn = std::bit_cast<uint32_t>(b);

	if (a >= 0 && b >= 0)
		return bPunn - aPunn;
	if (a < 0 && b < 0)
		return aPunn - bPunn;

	return (aPunn - 0x80000000) + bPunn;
}

std::vector<float> FloatRefine(const std::vector<double>& xs_, const std::vector<double>& dCoeffs)
{
	// === Parameters ===
	static constexpr size_t NumSamples = 1'000'000;
	static constexpr size_t NumIter = 1'000;
	// ==================

	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<float> branch{ 0.0f, 1.0f };
	std::uniform_real_distribution<float> dist{ (float)min, (float)max };

	// Round coefficients
	std::vector<float> bestCoeffs;
	for (double c : dCoeffs)
		bestCoeffs.push_back((float)c);

	uint32_t bestError = 10'000;
	for (size_t i = 0; i < NumIter; i++)
	{
		std::vector<float> newCoeffs{ bestCoeffs };
		for (float& coeff : newCoeffs)
		{
			if (branch(gen) < 0.1f)
				coeff = std::nextafter(coeff, (gen() & 1) ? -INFINITY : INFINITY);
		}

		// Calculate max error
		uint32_t maxError = 0;
		for (size_t i = 0; i < NumSamples; i++)
		{
			float x = dist(gen);
			float y = fFunc(x);

			float approx = EvaluateRational(x, newCoeffs);
			uint32_t err = ULPDistance(approx, y);

			if (err > maxError)
				maxError = err;
		}

		if (maxError < bestError)
		{
			bestError = maxError;
			bestCoeffs = newCoeffs;
			std::cout << "Error: " << maxError << '\n';
		}
	}

	return bestCoeffs;
}

volatile bool keepRunning = true;
void StopThread()
{
	std::cin.get();
	keepRunning = false;
}

int main()
{
	// Prepare data
	std::vector<double> xs, ys;
	for (double x = min; x < max; x += step)
	{
		xs.push_back(x);
		ys.push_back(Func(x));
	}

	// Initial parameters
	std::vector<double> ps{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1  };
	std::vector<double> scales{
		0,
		152.91909281655765,
		1024.8955781203358,
		2455.09202423243,
		2528.414754708975,
		1066.824695457973,
		148.54000057530493,
		4.0539143590175835,
		152.91909351081017,
		1177.814723036858,
		3403.528339686068,
		4573.00254043458,
		2878.9137364837616,
		761.1479653139235,
		65.5632478027466
	};

	if (ps.size() != scales.size() || ps.size() != numCoeffs)
		throw;

	// Start thread
	std::thread stopThread{ StopThread };

	// Gradient descent
	for (size_t iter = 0; keepRunning; iter++)
	{
		// Calculate error and derivatives
		auto [initialError, derivs] = GetMaxError(xs, ys, ps, scales);

		if (iter % 1000 == 0)
			std::cout << std::format("Iter {} error: {}\n", iter, log10(initialError));

		// Take steps
		for (size_t pi = 0; pi < numCoeffs; pi++)
			ps[pi] -= derivs[pi] * 1e-14;
	}
	stopThread.join();

	// Get coefficients
	std::vector<double> coeffs = GetCoefficients(ps, scales);

	// Print numerator coefficients
	for (size_t i = 0; i < pOrder + 1; i++)
		std::cout << std::format("{}\n", coeffs[i]);

	std::cout << "\n\n";

	// Print denominator coefficients
	for (size_t i = pOrder + 1; i < numCoeffs; i++)
		std::cout << std::format("{}\n", coeffs[i]);
	std::cout << "1\n";

#if 0
	std::cout << "=== Refining ===\n";
	auto fCoeffs = FloatRefine(xs, coeffs);

	// Print numerator coefficients
	for (size_t i = 0; i < pOrder + 1; i++)
		std::cout << std::format("{}\n", fCoeffs[i]);

	std::cout << "\n\n";

	// Print denominator coefficients
	for (size_t i = pOrder + 1; i < numCoeffs; i++)
		std::cout << std::format("{}\n", fCoeffs[i]);
	std::cout << "1\n";
#endif

	std::cin.get();
}