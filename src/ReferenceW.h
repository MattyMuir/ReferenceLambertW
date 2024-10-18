#pragma once
#include <mpfr.h>

#include "Interval.h"

#define TRACK_BISECTIONS 0

class ReferenceW
{
private:
	enum class Sign { Negative, Positive, Inconclusive };

public:
	ReferenceW();
	~ReferenceW();

	Interval W0(double x);
	Interval Wm1(double x);

#if TRACK_BISECTIONS
	void LogBisectionStats() const;
#endif

private:
	mpfr_t xMpfr, low, high, temp0, temp1;
	mpfr_t m, yLowP0, yHighP0, yLowP1, yHighP1;

#if TRACK_BISECTIONS
	size_t numBisections = 0;
	size_t numInconclusive = 0;
#endif

	Sign GetMidpointSign(double x, mpfr_t midpoint, bool useHighPrec);
	Interval Bisection(double x, mpfr_t low, mpfr_t high, bool increasing);
};