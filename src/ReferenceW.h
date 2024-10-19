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
	mpfr_t expUp, expDown, wexpDown, wexpUp, numerator0, numerator1, wplus2, denominator1, frac1, denominator0, newW;

#if TRACK_BISECTIONS
	size_t numBisections = 0;
	size_t numInconclusive = 0;
#endif

	Sign GetMidpointSign(double x, mpfr_t midpoint, bool useHighPrec);
	Interval Bisection(double x, mpfr_t low, mpfr_t high, bool increasing);
	void HalleyW0(mpfr_t result, mpfr_t x, mpfr_t w, bool isUpper);
	static void LogUpDown(mpfr_t down, mpfr_t up, mpfr_t x);
	static void ExpUpDown(mpfr_t down, mpfr_t up, mpfr_t x);
};