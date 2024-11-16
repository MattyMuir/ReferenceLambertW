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

	Sign GetMidpointSign(float x, float midpoint, bool useHighPrec);

private:
	mpfr_t m, yLowP0, yHighP0, yLowP1, yHighP1;

	Intervalf Bisection(float x, float low, float high, bool increasing);
};