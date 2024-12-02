#pragma once
#include <utility>

#include <mpfr.h>

#include "Interval.h"
#include "Sign.h"

class ReferenceW
{
public:
	ReferenceW();
	~ReferenceW();

	Interval W0(double x);
	Interval Wm1(double x);

#if REFERENCEW_STATS
	double GetHighPrecRate() const;
	size_t GetMaxBisections() const;
	double GetAvgBisections() const;
#endif

private:
	mpfr_t m, yLowP0, yHighP0, yLowP1, yHighP1;

#if REFERENCEW_STATS
	size_t numEvals = 0, numHighPrec = 0, maxBisections = 0, totalBisections = 0;
#endif

	std::pair<double, double> Wm1Bracket(double x);
	Sign GetMidpointSign(double x, double midpoint, bool useHighPrec);
	Interval Bisection(double x, double low, double high, bool increasing);
};