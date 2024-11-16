#pragma once
#include <cfenv>

#include "rndutil.h"

template <typename Ty>
Ty HalleyW0(Ty x, Ty w, bool isUpper)
{
	/*
	result			= w - numerator0 / denominator0
	numerator0		= w e^w - x
	denominator0	= e^w (w + 1) - numerator1 / denominator1
	numerator1		= (w + 2)(w e^w - x)
	denominator1	= 2w + 2

	=== Upper Bound Case ===
	- numerator0	is POSITIVE and needs to be rounded DOWN
	- denominator0	is POSITIVE and needs to be rounded UP
	- numerator1	is POSITIVE and needs to be rounded DOWN
	- denominator1	is POSITIVE and needs to be rounded UP

	=== Lower Bound Case ===
	- numerator0	is NEGATIVE and needs to be rounded UP
	- denominator0	is POSITIVE and needs to be rounded UP
	- numerator1	is NEGATIVE and needs to be rounded DOWN
	- denominator1	is POSITIVE and needs to be rounded DOWN
	*/

	int rnd;

	// expUp and expDown
	auto [expDown, expUp] = ExpUpDown(w);

	// wexpDown
	Ty exp0 = (w > 0) ? expDown : expUp;
	Ty wexpDown = mul(w, exp0, FE_DOWNWARD);

	// wexpUp
	Ty exp1 = (w > 0) ? expUp : expDown;
	Ty wexpUp = mul(w, exp1, FE_UPWARD);

	// numerator0
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	Ty numerator0 = sub(isUpper ? wexpDown : wexpUp, x, rnd);

	// numerator1
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	Ty wplus2 = add(w, 2, rnd);
	Ty numerator1 = sub(wexpDown, x, FE_DOWNWARD);
	numerator1 = mul(numerator1, wplus2, FE_DOWNWARD);

	// denominator1
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	Ty denominator1 = add(w, w, rnd);
	denominator1 = add(denominator1, 2, rnd);

	// frac1
	Ty frac1 = div(numerator1, denominator1, FE_DOWNWARD);

	// denominator0
	Ty denominator0 = add(w, 1, FE_UPWARD);
	denominator0 = mul(denominator0, expUp, FE_UPWARD);
	denominator0 = sub(denominator0, frac1, FE_UPWARD);

	// newW
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	Ty newW = div(numerator0, denominator0, rnd);
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	newW = sub(w, newW, rnd);

	return newW;
}

template <typename Ty>
Ty HalleyWm1(Ty x, Ty w, bool isUpper)
{
	/*
	result			= w - numerator0 / denominator0
	numerator0		= w e^w - x
	denominator0	= e^w (w + 1) - numerator1 / denominator1
	numerator1		= (w + 2)(w e^w - x)
	denominator1	= 2w + 2

	=== Upper Bound Case ===
	- numerator0	is NEGATIVE				and needs to be rounded UP
	- denominator0	is NEGATIVE				and needs to be rounded DOWN
	- numerator1	is NEGATIVE when w > -2 and needs to be rounded DOWN
	- denominator1	is NEGATIVE				and needs to be rounded UP when w > -2

	=== Lower Bound Case ===
	- numerator0	is POSITIVE				and needs to be rounded DOWN
	- denominator0	is NEGATIVE				and needs to be rounded DOWN
	- numerator1	is POSITIVE when w > -2 and needs to be rounded DOWN
	- denominator1	is NEGATIVE				and needs to be rounded DOWN when w > -2
	*/

	int rnd;

	// expUp and expDown
	auto [expDown, expUp] = ExpUpDown(w);

	// wexpDown
	Ty wexpDown = mul(w, expUp, FE_DOWNWARD);

	// wexpUp
	Ty wexpUp = mul(w, expDown, FE_UPWARD);

	// numerator0
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	Ty numerator0 = sub(isUpper ? wexpUp : wexpDown, x, rnd);

	// numerator1
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	Ty wplus2 = add(w, 2, rnd);
	rnd = (wplus2 > 0) ? FE_DOWNWARD : FE_UPWARD;
	Ty numerator1 = sub((wplus2 > 0) ? wexpDown : wexpUp, x, rnd);
	numerator1 = mul(numerator1, wplus2, FE_DOWNWARD);

	// denominator1
	rnd = (isUpper == (wplus2 > 0)) ? FE_UPWARD : FE_DOWNWARD;
	Ty denominator1 = add(w, w, rnd);
	denominator1 = add(denominator1, 2, rnd);

	// frac1
	Ty frac1 = div(numerator1, denominator1, FE_UPWARD);

	// denominator0
	Ty denominator0 = add(w, 1, FE_DOWNWARD);
	denominator0 = mul(denominator0, expUp, FE_DOWNWARD);
	denominator0 = sub(denominator0, frac1, FE_DOWNWARD);

	// newW
	rnd = isUpper ? FE_DOWNWARD : FE_UPWARD;
	Ty newW = div(numerator0, denominator0, rnd);
	rnd = isUpper ? FE_UPWARD : FE_DOWNWARD;
	newW = sub(w, newW, rnd);

	return newW;
}