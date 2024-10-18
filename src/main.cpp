#include <iostream>
#include <iomanip>

#include "ReferenceLambertW.h"

int main()
{
	double x = 0;
	Interval result = ReferenceW0(x);
	std::cout << std::setprecision(17);
	std::cout << "W(" << x << ") = [" << result.inf << ", " << result.sup << "]\n";
}