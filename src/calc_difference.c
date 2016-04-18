#include <math.h>
#include <stdlib.h>
#include "misc_functions.h"
/* 
	Implementation by: Tallulah Andrews
	Date: 15 April 2016
*/

void calc_difference (int* expr_mat, int* nc, int* ng, double* difference) {
	long tis[*nc]; 
	long sj[*ng];
	long total;
	
	int i,j;
	for (i = 0; i < *nc; i++) {
		for (j = 0; j < *ng; j++) {
			int val = expr_mat[convert_2D_indices_to_1D(j, i, ng, nc)];
			tis[i] += val;
			sj[j] += val;
			total += val;
		}
	}
	for (i = 0; i < *nc; i++) {
		for (j = 0; j < *ng; j++) {
			double lambda = tis[i]*sj[j]/((double)total);
			int coord = convert_2D_indices_to_1D(j, i, ng, nc);
			int val = expr_mat[coord];
			difference[coord] = (double)val - lambda;
		}
	}	
}
