#pragma once
#include <random>
#include <numeric>

// Reciprocal distribution between two positive, finite endpoints
template <typename Ty>
class ReciprocalDistribution
{
public:
	ReciprocalDistribution(Ty min_, Ty max_)
		: min(min_), max(max_), dist(0, 1)
	{
		int minExp, maxExp;
		Ty minMan = std::frexp(min, &minExp);
		Ty maxMan = std::frexp(max, &maxExp);

		// Compute logWidth
		if (maxExp - minExp > 10)
			logWidth = std::log(max) - std::log(min);
		else
			logWidth = std::log1p((max - min) / min);

		// Compute geometricMean
		Ty meanMan = sqrt(minMan * maxMan);
		int expSum = minExp + maxExp;
		geometricMean = std::ldexp(meanMan, expSum / 2);
		if (expSum % 2)
			geometricMean *= (Ty)1.4142135623730950488;

		// Check if mean is needed
		Ty expMidpoint = exp(logWidth / 2);
		useMean = !std::isfinite(min * expMidpoint) || !std::isfinite(max / expMidpoint);
	}

	Ty Min() const
	{
		return min;
	}
	Ty Max() const
	{
		return max;
	}

	template <typename Engine>
	Ty operator()(Engine& eng)
	{
		Ty val = dist(eng);

		if (useMean && val > 0.25 && val < 0.75)
			return geometricMean * std::exp((val - (Ty)0.5) * logWidth);

		if (val < 0.5)
			return min * std::exp(val * logWidth);
		return max / std::exp((1 - val) * logWidth);
	}

protected:
	Ty min, max, logWidth, geometricMean;
	bool useMean;
	std::uniform_real_distribution<Ty> dist;

	ReciprocalDistribution() = default;

	template <typename Ty_>
	friend class ReciprocalDistributionEx;
};