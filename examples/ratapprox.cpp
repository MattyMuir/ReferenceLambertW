#include <cassert>

#include <iostream>
#include <vector>
#include <format>
#include <random>

#include <ReferenceLambertW.h>

constexpr size_t pOrder = 8;
constexpr size_t qOrder = 7;
constexpr size_t numCoeffs = pOrder + qOrder + 1;
constexpr double yOffset = 0.567143290409783872999968662210;

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
	return W0(exp(x));
}

float fFunc(float x)
{
	static ReferenceWf evaluator;
	return evaluator.Wm1(-exp(-0.5f * (x * x + 2))).inf;
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
	// Round xs
	std::vector<float> xs;
	for (double x : xs_)
		xs.push_back((float)x);

	// Calculate ys
	std::vector<float> ys;
	for (float x : xs)
		ys.push_back((float)Func(x));

	// Round coefficients
	std::vector<float> bestCoeffs;
	for (double c : dCoeffs)
		bestCoeffs.push_back((float)c);

	uint32_t bestError = 10'000;
	static std::mt19937_64 gen{ std::random_device{}() };
	static std::uniform_real_distribution<float> branch{ 0.0f, 1.0f };
	static std::normal_distribution<float> dist{ 1, 1e-5 };
	for (size_t i = 0; i < 1'000'000; i++)
	{
		std::vector<float> newCoeffs{ bestCoeffs };
		for (float& coeff : newCoeffs)
		{
			if (branch(gen) < 0.1f)
				coeff = coeff * dist(gen);
		}

		// Calculate max error
		uint32_t maxError = 0;
		for (size_t i = 0; i < xs.size(); i++)
		{
			float x = xs[i];
			float y = ys[i];

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

int main()
{
	// Prepare data
	std::vector<double> xs, ys;
	for (double x = 0; x < 79.52; x += 0.05)
	{
		xs.push_back(x);
		ys.push_back(Func(x));
	}

	// Initial parameters
	std::vector<double> ps{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	std::vector<double> scales{
		0,
		11087152.7748,
		8486181.12685,
		3416914.71516,
		813083.430535,
		115397.465193,
		9766.42319992,
		230.594326831,
		0.998694488778,
		30636450.9383,
		17210757.1717,
		6054509.95539,
		1212699.80199,
		155314.797354,
		10977.9615602,
		237.729006177
	};

	if (ps.size() != scales.size() || ps.size() != numCoeffs)
		throw;

	// Gradient descent
	for (size_t iter = 0;; iter++)
	{
		// Calculate error and derivatives
		auto [initialError, derivs] = GetMaxError(xs, ys, ps, scales);

		// Stop condition
		if (log10(initialError) < -9)
			break;

		if (iter % 1000 == 0)
			std::cout << std::format("Iter {} error: {}\n", iter, log10(initialError));

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
}