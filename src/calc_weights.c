/* Copyright (c) 2016 Genome Research Ltd .
Author : Tallulah Andrews <tallulandrews@gmail.com>
This file is part of PoissonUMIs.

PoissonUMIs is free software : you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program . If not , see <http://www.gnu.org/licenses/>. */

#include <math.h>
#include "misc_functions.h"

double poisson_cdf_weight (double val, double lambda) {
	// val = obs*alpha, lambda = sj*ti*totalreads*alpha
	double diff = fabs(val-lambda);
	int val_ceil = (int) ceil( (double)(lambda-diff) );
	int val_floor = (int) floor( (double)(lambda+diff) );
	double phigh = pPoisson( val_floor, lambda );
	double plow = pPoisson( val_ceil, lambda );
	double total = phigh+plow;
	if (total > 1.0) { //might be possible since pPoisson is not exact, or if diff == 0 and lambda is an integer.
		total = 1.0;
	}
	return(total);
}

void calc_weights (double* val_mat, double* lambda_mat, double* weight_mat, int* ncol, int* nrow) {
	int row, col;
	for (row = 0; row < *nrow; row++) {
		for (col = 0; col < *ncol; col++) {
			int coord = convert_2D_indices_to_1D(row, col, nrow, ncol);
			weight_mat[coord] = 1.0 - poisson_cdf_weight( val_mat[coord], lambda_mat[coord] );
		}
	}
}
