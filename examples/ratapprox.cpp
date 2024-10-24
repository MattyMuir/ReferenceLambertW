#include <iostream>
#include <vector>
#include <format>

#include "ReferenceW.h"

constexpr size_t pOrder = 5;
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
	for (double x = 0.5; x < 709; x += 0.1)
	{
		xs.push_back(x);
		ys.push_back(evaluator.Wm1(-exp(-x - 1)).inf);
	}

	// Initial parameters
	std::vector<double> ps{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<double> scales{
		-168122.526858,
		-598634.985232,
		-304540.505252,
		-29876.9519565,
		-493.433883072,
		-1.00044283758,
		129082.609542,
		193651.02343,
		26939.8181948,
		484.946505393
	};

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
			static constexpr double h = 1e-8;
			ps[pi] += h;
			double newError = GetMaxError(xs, ys, ps, scales);
			derivs[pi] = (newError - initialError) / h;

			ps[pi] = originalParam;
		}

		// Take steps
		for (size_t pi = 0; pi < numCoeffs; pi++)
			ps[pi] -= derivs[pi] * 1e-5;
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