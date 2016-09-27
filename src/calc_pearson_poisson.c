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
#include <stdlib.h>
#include "misc_functions.h"
/* 
	Implementation by: Tallulah Andrews
	Date: 2 May 2016
*/

void pearson_cor_poisson (int* expr_mat, int* nc, int* ng, double* lambdas, double* corrs) {
	// All pairs of columns
	int col1, col2;
	int row;
	for (col1 = 0; col1 < (*nc-1); col1++) {
		for(col2 = col1+1; col2 < *nc; col2++) {
			double var1 = 0.0; 
			double var2 = 0.0;
			double covar = 0.0;
			for (row = 0; row < *ng; row++) {
				int coords_x = convert_2D_indices_to_1D(row, col1, ng, nc);
				int coords_y = convert_2D_indices_to_1D(row, col2, ng, nc);
				covar += (expr_mat[coords_x] - lambdas[coords_x])*(expr_mat[coords_y] - lambdas[coords_y]);
				var1 += (expr_mat[coords_x]-lambdas[coords_x])*(expr_mat[coords_x]-lambdas[coords_x]);
				var2 += (expr_mat[coords_y]-lambdas[coords_y])*(expr_mat[coords_y]-lambdas[coords_y]);
			}
			int coord_out1 = convert_2D_indices_to_1D(col1, col2, nc, nc);
			int coord_out2 = convert_2D_indices_to_1D(col2, col1, nc, nc);
			double cor = covar/(sqrt(var1)*sqrt(var2));
			corrs[coord_out1] = cor;
			corrs[coord_out2] = cor;
		}
	}
}
