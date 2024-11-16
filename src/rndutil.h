#pragma once
#include <utility>

#include <mpfr.h>

// float
float add(float x, float y, int rnd);
float sub(float x, float y, int rnd);
float mul(float x, float y, int rnd);
float div(float x, float y, int rnd);
float sqrt(float x, int rnd);

std::pair<float, float> ExpUpDown(float x);
std::pair<float, float> LogUpDown(float x);

// double
double add(double x, double y, int rnd);
double sub(double x, double y, int rnd);
double mul(double x, double y, int rnd);
double div(double x, double y, int rnd);
double sqrt(double x, int rnd);

std::pair<double, double> ExpUpDown(double x);
std::pair<double, double> LogUpDown(double x);

// mpfr
void ExpUpDown(mpfr_t down, mpfr_t up, mpfr_t x);