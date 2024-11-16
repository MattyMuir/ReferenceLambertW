#pragma once
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

	Sign GetMidpointSign(double x, double midpoint, bool useHighPrec);

private:
	mpfr_t m, yLowP0, yHighP0, yLowP1, yHighP1;

	Interval Bisection(double x, double low, double high, bool increasing);
};