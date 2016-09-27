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
#define ITMAX 100 //max iterations
#define EPS 1.0e-3 //accuracy
#define FPMIN 1.0e-30 //Number near smallest possilbe doubleing-point number.

int convert_2D_indices_to_1D (int i, int j, int* nrow, int* ncol) {
	return(j* (*nrow) + i);
}

/* Returns ln(gamma(xx)) for xx > 0  with error < 2*10^-10 
From: Numerical Recipies in C.
*/
double gammln(double xx) {
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146, -86.50532032941677,
			      24.01409824083091, -1.231739572450155,
			      0.120865097386617e-2, -0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005+ser/x);
}

double dPoisson (int x, double lambda) {
	if (x == 0) {
		return(exp(-lambda));
	}
	double lnp = x*log(lambda) - lambda - (gammln(x+1.0));
	return(exp(lnp));
}

/* Poisson cdf from Numerical Recipes in C */
double pPoisson (int k, double lambda) {
	double gamser, gammcf, gln;

//	if (lambda < 0.0 || k < 0) nrerror("Invalid arguments in routine gammq");
	if (k == 0) {
		return(dPoisson(k, lambda));
	}
	if (k < 25) {
		double sum = 0.0;
		int i;
		for (i = 0; i <=k; i++) {
			sum+=dPoisson(i, lambda);
		}
		return(sum);
	}
	if (lambda < (k + 1.0)) {
		gser(&gamser, (double) k, lambda, &gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf, (double) k, lambda, &gln);
		return gammcf;
	}
}

void gser(double* gamser, double a, double x, double* gln) {
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
//		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1; n<=ITMAX; n++) {
			++ap;
			del += x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		//nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}

void gcf(double* gammcf, double a, double x, double* gln) {
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h += del;
		if (fabs(del-1.0) < EPS) break;
	}
//	if (i > ITMAX) nerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

