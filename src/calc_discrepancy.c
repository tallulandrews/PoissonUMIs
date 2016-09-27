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
	Date: 15 April 2016
*/

void calc_discrepancy (int* expr_mat, int* nc, int* ng, double* discrepancy) {
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
			if (val > lambda) {
				discrepancy[coord] = pPoisson(val-1,lambda);
			} else {
				discrepancy[coord] = -1*(1-pPoisson(val,lambda));
			}
		}
	}	
}
