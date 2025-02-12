#include <iostream>
#include <vector>
#include <format>
#include <random>
#include <thread>

#include <ReferenceLambertW.h>

// === Parameters ===
constexpr size_t pOrder = 3;
constexpr size_t qOrder = 3;
constexpr size_t numCoeffs = pOrder + qOrder + 1;
constexpr double yOffset = 0;

constexpr double min = -0.28;
constexpr double max = 7.34;
constexpr double step = 0.0051;
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
	//return Wm1(-exp(-x * x - 1));
	//return W0(exp(x));
	//return W0(exp(x) - 1);
	//return W0(x);
	return (W0(x) - x) / (x * x);
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
	std::vector<double> ps{ 1, 1, 1, 1, 1, 1, 1 };
	std::vector<double> scales{
		-0.811183167505,
		-2.09065201806,
		-0.919464181645,
		-0.00253967039523,
		0.811104259067,
		3.30625740088,
		3.71791903293
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
			ps[pi] -= derivs[pi] * 1e-8;
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
}