#include <stdlib.h>
#include "misc_functions.h"

double* calc_lambdas (int* expr_mat, int* nc, int* ng) {
	double* lambdas = (double*)calloc((*nc) * (*ng),  sizeof(double));
	if (NULL == lambdas) { //possible error
		return(lambdas);
	}

	long tis[*nc]; //total reads per cell
	long sj[*ng]; //total reads per gene
	long total; //total reads in matrix
	
	int i,j;
	for (i = 0; i < *nc; i++) {
		for (j =0; j < *ng; j++) {
			int val = expr_mat[convert_2D_indices_to_1D(j,i,ng,nc)];
			tis[i] += val;
			sj[j] += val;
			total += val;
		}
	}
	for (i = 0; i < *nc; i++) {
		for (j =0; j < *ng; j++) {
			double lambda = tis[i]*sj[j]/((double)total);
			lambdas[convert_2D_indices_to_1D(j,i,ng,nc)] = lambda;
		}
	}
	return(lambdas);
}	
