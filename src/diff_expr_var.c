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
#include <stdio.h>
#include <stdlib.h>
#include "misc_functions.h"
#include "poilog_package_C_code.h"

double calc_gene_log_likelihood (double* val_mat, double* lambda_mat, int* ncol, int* nrow, int gene_row) {
	int col;
	double likelihood = 1.0;
	for (col = 0; col < *ncol; col++) {
		int coord = convert_2D_indices_to_1D(gene_row, col, nrow, ncol);
		double rounded_val = ceil(val_mat[coord]);
		double prob = dPoisson( (int)rounded_val, lambda_mat[coord] );
		if (prob < 0.0001) {prob = 0.0001;}
		likelihood = likelihood + log(prob);
	}
	return(likelihood);
}

double chi_sq_test (double chi_sq, int df) {
	return( pPoisson(df/2, chi_sq/2.0) ); //Is integer division a problem here?
}

/* output should be a nrow x nlabels+1 matrix*/
/* labels should be 0-X */
void diff_expression (double* y, double* lambdas, int* nrow, int* ncol, int* labels, int* nlabels, double* out) {
	int row1, column; 
//	double group_lambdas[*nrow * *ncol];

	/* Calculate lambdas where each group has it's own mean expression for 
	the gene. 
	This is calculated by scaling the existing global lambda by a factor
	which represents the mean expression in the group relative to the
	global mean expression.*/
	double* group_lambdas = (double*)malloc(*nrow * *ncol * sizeof(double));
	int* ncol_out = (int*)malloc(sizeof(int));
	*ncol_out = *nlabels+1;
	for (row1 = 0; row1 < *nrow; row1++) {

		/* Calculate expression & size of each group */
		long group_sums[*nlabels]; //total molecules per group
		long group_lambda_sums[*nlabels]; //total molecules per group
		int group_size[*nlabels]; //number of cells per group
		long row_sum = 0; //total molecules for this gene across all cells
		int row_size = 0;

		// initialize
		int tmp;
		for (tmp = 0; tmp < *nlabels; tmp++) {
			group_size[tmp] = 0;
			group_sums[tmp] = 0;
			group_lambda_sums[tmp] = 0;
		}

		for (column = 0; column < *ncol; column++) { //across cells
			int index = convert_2D_indices_to_1D(row1, column, nrow, ncol);
			double val = y[index];
			val = ceil(val);
			row_sum = row_sum + val;
			row_size = row_size + 1;
			int group_id = labels[column];
			group_sums[group_id] = group_sums[group_id] + val;
			group_lambda_sums[group_id] = group_lambda_sums[group_id] + lambdas[index];
			group_size[group_id] = group_size[group_id] + 1;
		}
		/* Scale each lambda by the relevant group factor */
		for (column = 0; column < *ncol; column++) { 
			int index = convert_2D_indices_to_1D(row1, column, nrow, ncol);
			int group_id = labels[column];
			double group_factor = 1.0;
			if (group_sums[group_id] > 0) {
				group_factor = ((double)group_sums[group_id] / group_size[group_id] ) / ( (double)row_sum / (row_size)); //relative mean expression
			} else {
				if (group_size[group_id] > 0) {
					group_factor = 1.0/ (double)group_size[group_id];
				}
			}

			group_lambdas[index] = lambdas[index] * group_factor;
			out[convert_2D_indices_to_1D(row1, group_id, nrow, ncol_out)] = ( (double)group_sums[group_id] - (double) group_lambda_sums[group_id] )/group_size[group_id]; //store mean expression for each group in output
		}
		double likelihood_null = calc_gene_log_likelihood(y, lambdas, ncol, nrow, row1);
		double likelihood_alt = calc_gene_log_likelihood(y, group_lambdas, ncol, nrow, row1);
		double log_ratio = likelihood_alt-likelihood_null;
		if (log_ratio < 0) {log_ratio=0;}
		double stat = 2.0*log_ratio;
		double pval = chi_sq_test(stat, (*nlabels)-1);
		out[convert_2D_indices_to_1D(row1, (*nlabels), nrow, ncol_out)] = pval;
	}
}

void distance_wgt (double* y, double* w, int* nrow, int* ncol, double* exponent, double* out) {
	int row1,row2, column;
	for (row1 = 0; row1 < *nrow; row1++) {
		for (row2 = row1; row2 < *nrow; row2++) {
			double dist = 0.0;
			double sum_w = 0.0;
			if (row2 == row1) {
				out[convert_2D_indices_to_1D(row1, row2, nrow, nrow)] = 0.0;
				continue;
			}
			for ( column = 0; column < *ncol; column++) {
				int coords1 = convert_2D_indices_to_1D(row1, column, nrow, ncol);
				int coords2 = convert_2D_indices_to_1D(row2, column, nrow, ncol);
				double difference = fabs(y[coords1]-y[coords2]);
//				double weight = w[coords1]*w[coords2];
				double weight = w[coords1]+w[coords2];
				sum_w = sum_w + weight;
				
				if (*exponent == 2.0) {
					dist = dist + weight*difference*difference;
				} else if (*exponent == 1.0) {
					dist = dist + weight*difference;
				} else {
					dist = dist + weight*pow(difference, *exponent);
				}
			}
			
			dist = dist/sum_w;
			if (*exponent == 2.0) {
                                dist = sqrt(dist);
                        } else if (*exponent == 1.0) {
                                //dist = dist;
                        } else {
                                dist = pow(dist, 1.0/(*exponent));
                        }
			out[convert_2D_indices_to_1D(row1, row2, nrow, nrow)] = dist;
			out[convert_2D_indices_to_1D(row2, row1, nrow, nrow)] = dist;
		}
	}
}
