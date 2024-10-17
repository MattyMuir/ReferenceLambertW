#include <iostream>
#include <iomanip>

#include "ReferenceLambertW.h"

int main()
{
	std::cout << std::setprecision(20);

	freopen("out.csv", "w", stdout);
	for (double x = log(20); x < 709; x += 5)
	{
		Interval res = ReferenceW0(exp(x));

		std::cout << x << ',' << res.inf << '\n';
	}
}