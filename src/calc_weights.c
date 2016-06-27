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
