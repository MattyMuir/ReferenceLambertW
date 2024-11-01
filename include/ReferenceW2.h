#pragma once
#include <mpfr.h>

#include "../src/Interval.h"

#define TRACK_BISECTIONS 0

class ReferenceW2
{
public:
	enum class Sign { Negative, Positive, Inconclusive };

public:
	ReferenceW2();
	~ReferenceW2();

	Interval W0(double x);
	Interval Wm1(double x);

#if TRACK_BISECTIONS
	void LogBisectionStats() const;
#endif

	Sign GetMidpointSign(double x, double m, bool useHighPrec);

private:
	mpfr_t m, yLowP0, yHighP0, yLowP1, yHighP1;

#if TRACK_BISECTIONS
	size_t numBisections = 0;
	size_t numInconclusive = 0;
#endif

	Interval Bisection(double x, double low, double high, bool increasing);
	double HalleyW0(double x, double w, bool isUpper);
	double HalleyWm1(double x, double w, bool isUpper);
};