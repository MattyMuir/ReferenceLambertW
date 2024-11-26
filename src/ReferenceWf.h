#pragma once
#include <mpfr.h>

#include "Interval.h"
#include "Sign.h"

class ReferenceWf
{
public:
	ReferenceWf();
	~ReferenceWf();

	Intervalf W0(float x);
	Intervalf Wm1(float x);

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

	Sign GetMidpointSign(float x, float midpoint, bool useHighPrec);
	Intervalf Bisection(float x, float low, float high, bool increasing);
};