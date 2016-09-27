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

int convert_2D_indices_to_1D (int i, int j, int* nrow, int* ncol);

double gammln(double xx);

double dPoisson(int x, double lambda);

double pPoisson(int k, double lambda);

void gser(double* gamser, double a, double x, double* gln);

void gcf(double* gammcf, double a, double x, double* gln);

//void nrerror(char error_text[]);

