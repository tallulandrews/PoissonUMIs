#include <math.h>
#include <stdlib.h>
#include "misc_functions.h"
/* 
	Implementation by: Tallulah Andrews
	Date: 5 Oct 2015
	Untested
*/

void covariance_weighted (double* x, double* wx, double* y, double* wy, int* n, double* covar) {
	covar[0] = (double) *n;
	double sum_x = 0.0;
	double sum_wx = 0.0;
	double sum_y = 0.0;
	double sum_wy = 0.0;
	int i;
	for (i = 0; i < *n; i++) {
		sum_x = sum_x + wx[i]*x[i];
		sum_wx = sum_wx + wx[i];
		sum_y = sum_y + wy[i]*y[i];
		sum_wy = sum_wy + wy[i];
	}
	double sum_w = 0.0;
	double sum_cov = 0.0;
	for (i = 0; i < *n; i++) {
		double w = wx[i]*wy[i];
		sum_w = sum_w + w;
		sum_cov = sum_cov + w*(x[i]-sum_x/sum_wx)*(y[i]-sum_y/sum_wy);
	}
	covar[0] = 1.0/(sum_w-1.0)*sum_cov;
}



void cov_matrix_weighted (double* m, double* w, int* nrow, int* ncol, double* out) {
	int row1,row2;
	for (row1 = 0; row1 < *nrow; row1++) {
		for (row2 = row1; row2 < *nrow; row2++) {
			double covar = 0.0;
			int coord_row1 = convert_2D_indices_to_1D(row1,0,nrow,ncol);
			int coord_row2 = convert_2D_indices_to_1D(row2,0,nrow,ncol);
			covariance_weighted (&m[coord_row1], &w[coord_row1], &m[coord_row2], &w[coord_row2], ncol, &covar);
			out[convert_2D_indices_to_1D(row1,row2,nrow,nrow)] = covar;
			out[convert_2D_indices_to_1D(row2,row1,nrow,nrow)] = covar;
		}
	}	
}

void cor_matrix_weighted (double* m, double* w, int* nrow, int* ncol, double* out) {
	int row1,row2;
	for (row1 = 0; row1 < *nrow; row1++) {
		double var1 = 0.0;
		int coord_row1 = convert_2D_indices_to_1D(row1,0,nrow,ncol);
		covariance_weighted (&m[coord_row1], &w[coord_row1], &m[coord_row1], &w[coord_row1], ncol, &var1);
		for (row2 = row1; row2 < *nrow; row2++) {
			double covar = 0.0;
			double var2 = 0.0;
			int coord_row2 = convert_2D_indices_to_1D(row2,0,nrow,ncol);
			covariance_weighted (&m[coord_row1], &w[coord_row1], &m[coord_row2], &w[coord_row2], ncol, &covar);
			covariance_weighted (&m[coord_row2], &w[coord_row2], &m[coord_row2], &w[coord_row2], ncol, &var2);
			out[convert_2D_indices_to_1D(row1,row2,nrow,nrow)] = covar/(sqrt(var1)*sqrt(var2));
			out[convert_2D_indices_to_1D(row2,row1,nrow,nrow)] = covar/(sqrt(var1)*sqrt(var2));
		}
	}	
}
