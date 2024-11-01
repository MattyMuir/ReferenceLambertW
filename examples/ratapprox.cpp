#include <cassert>

#include <iostream>
#include <vector>
#include <format>

#include "ReferenceW.h"

constexpr size_t pOrder = 4;
constexpr size_t qOrder = 4;
constexpr size_t numCoeffs = pOrder + qOrder + 1;

std::vector<double> GetCoefficients(const std::vector<double>& ps, const std::vector<double>& scales)
{
	std::vector<double> coeffs = ps;
	for (size_t i = 0; i < coeffs.size(); i++)
		coeffs[i] *= scales[i];

	return coeffs;
}

double GetMaxError(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& ps, const std::vector<double>& scales)
{
	// Calculate coefficients
	std::vector<double> coeffs = GetCoefficients(ps, scales);

	double maxError = 0;
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
		double approx = numer / denom;
		double error = abs((approx - y) / y);
		if (error > maxError)
			maxError = error;
	}

	return maxError;
}

int main()
{
	// Prepare data
	std::vector<double> xs, ys;
	ReferenceW evaluator;
	for (double x = -0.29; x < 7.34; x += 0.011)
	{
		xs.push_back(x);
		ys.push_back(evaluator.W0(x).inf);
	}

	// Initial parameters
	std::vector<double> ps{ 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<double> scales{
		0,
		30.6615230228,
		83.8519914551,
		46.1167994941,
		3.47179067097,
		30.6666937031,
		114.471438222,
		114.536515609,
		28.6775252301
	};

	if (ps.size() != scales.size() || ps.size() != numCoeffs)
		throw;

	// Gradient descent
	std::vector<double> derivs(numCoeffs);
	for (size_t iter = 0;; iter++)
	{
		// Calculate initial error
		double initialError = GetMaxError(xs, ys, ps, scales);

		// Stop condition
		if (log10(initialError) < -3.8)
			break;

		if (iter % 100 == 0)
			std::cout << std::format("Iter {} error: {}\n", iter, log10(initialError));

		// Calculate derivatives
		for (size_t pi = 0; pi < numCoeffs; pi++)
		{
			double originalParam = ps[pi];

			// Finite difference gradient approximation
			static constexpr double h = 1e-10;
			ps[pi] += h;
			double newError = GetMaxError(xs, ys, ps, scales);
			derivs[pi] = (newError - initialError) / h;

			ps[pi] = originalParam;
		}

		// Take steps
		for (size_t pi = 0; pi < numCoeffs; pi++)
			ps[pi] -= derivs[pi] * 1e-8;
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