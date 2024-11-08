#include <cassert>

#include <iostream>
#include <vector>
#include <format>

#include "ReferenceW.h"

constexpr size_t pOrder = 3;
constexpr size_t qOrder = 1;
constexpr size_t numCoeffs = pOrder + qOrder + 1;
constexpr double yOffset = -1;

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
		double y = ys[i];

		// Evaluate rational
		double numer = coeffs[pOrder];
		for (size_t i = 0; i < pOrder; i++)
			numer = numer * x + coeffs[pOrder - 1 - i];

		double denom = 1;
		for (size_t i = 0; i < qOrder; i++)
			denom = denom * x + coeffs[pOrder + qOrder - i];

		// Calculate error
		double approx = numer / denom + yOffset;
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
	return Wm1(-exp(-0.5 * x * x - 1));
}

int main()
{
	// Prepare data
	std::vector<double> xs, ys;
	for (double x = 0.1; x < 38; x += 0.1)
	{
		xs.push_back(x);
		ys.push_back(Func(x));
	}

	// Initial parameters
	std::vector<double> ps{ 1, 1, 1, 1, 1 };
	std::vector<double> scales{
		0,
		-5.5247349255,
		-2.83789677377,
		-0.499408296334,
		5.5247349255
	};

	if (ps.size() != scales.size() || ps.size() != numCoeffs)
		throw;

	// Gradient descent
	for (size_t iter = 0;; iter++)
	{
		// Calculate error and derivatives
		auto [initialError, derivs] = GetMaxError(xs, ys, ps, scales);

		// Stop condition
		if (log10(initialError) < -4.17)
			break;

		if (iter % 1000 == 0)
			std::cout << std::format("Iter {} error: {}\n", iter, log10(initialError));

		// Take steps
		for (size_t pi = 0; pi < numCoeffs; pi++)
			ps[pi] -= derivs[pi] * 1e-6;
	}

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