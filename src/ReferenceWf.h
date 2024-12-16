#pragma once
#include <utility>

#include <arb.h>

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
	arb_t xArb, mArb, yArb;

#if REFERENCEW_STATS
	size_t numEvals = 0, numHighPrec = 0, maxBisections = 0, totalBisections = 0;
#endif

	static std::pair<float, float> W0Bracket(float x);
	static std::pair<float, float> Wm1Bracket(float x);
	Sign GetMidpointSign(float x, float midpoint, bool useHighPrec);
	Intervalf Bisection(float x, float low, float high, bool increasing);
};