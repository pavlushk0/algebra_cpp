
#include <iostream>
#include <string>
#include <algorithm>
#include <tuple>

#include "algebra.h"

int main(int argc, char *argv[]) {
	mtrx3_t m, lm, um, mm;

	m = mtrx3_set(-10.0, 1.0, 5.0,
				   2.0, 7.0, 4.0,
			      -8.0, 4.0, 1.0);

	mtrx3_show(m);
	printf("\n");

	auto ret = mtrx3_lu(m);

	lm = std::get<0>(ret);
	mtrx3_show(lm);
	printf("\n");

	um = std::get<1>(ret);
	mtrx3_show(um);
	printf("\n");

	mm = mtrx3_mult(um, lm);

	mtrx3_show(mm);
	printf("\n");

	return 0;
}
