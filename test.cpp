
#include <iostream>
#include <string>
#include <algorithm>
#include <tuple>

#include "algebra.h"

int main(int argc, char *argv[]) {
	mtrx4_t m, lm, um, mm;

	m = mtrx4_t(-10.0f, 1.0f, 5.0f, 1,
				  2.0f, 7.0f, 4.0f, 1,
			     -8.0f, 4.0f, 1.0f, 1,
				 3, 4, 5, 6);

	mtrx_show<mtrx4_t, 4>(m);
	printf("\n");

	auto ret = mtrx_lu<mtrx4_t, 4>(m);

	lm = std::get<0>(ret);
	mtrx_show<mtrx4_t, 4>(lm);
	printf("\n");

	um = std::get<1>(ret);
	mtrx_show<mtrx4_t, 4>(um);
	printf("\n");

	mm = mtrx_mult<mtrx4_t, 4>(um, lm);

	mtrx_show<mtrx4_t, 4>(mm);
	printf("\n");

	return 0;
}
