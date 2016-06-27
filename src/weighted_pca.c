#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <R.h>
#include "misc_functions.h"

/* 
	Initialization procedure to minimize chance of sub-optimal convergence from: Gabriel, K. Ruben, and Zamir, S., (1979) Lower Rank Approximation of Matrices by Least Squares With Any Choice of Weights. Technometrics: 21(4) p489-498. 
	Generalization to multiple components from: Tamuz, O., Mazeh, T., and Zucker, S. (2013) Correcting systematic effects in a large set of photometric lightcurves. Mon. Not. R, Astron. Soc. 

	Implementation by: Tallulah Andrews
	Date: 28 Sept 2015
	Untested


*/

void initialize (double* y, double* w, int* Ng, int* Nc, double* a) {
	double close_enough = 0.00001; // close enough to 0 to count as a zero for sub-optimal convergence problems
	int exists_zero = 0;
	double max_Phi = -1;
	double max_Theta = -1; // weighted norm
	int maxPhi_column = -1; // column of max
	int maxTheta_column = -1; // column of max
	
	int i,j;
	for (i = 0; i < *Ng; i++) {
		for (j = 0; j < *Nc; j++) {
			if (w[convert_2D_indices_to_1D(i,j,Ng,Nc)] < close_enough) {
				exists_zero = 1;
				// Calculate Phi
				double sum1 = 0.0;
				double sum2 = 0.0;
				int e,g;
				for (e = 0; e < *Ng; e++) {
					if (e != i) { sum1 = sum1 + w[convert_2D_indices_to_1D(e,j,Ng,Nc)]*y[convert_2D_indices_to_1D(e,j,Ng,Nc)]*y[convert_2D_indices_to_1D(e,j,Ng,Nc)];}
				}
				for (g = 0; g < *Nc; g++) {
					if (g != j) { sum2 = sum2 + w[convert_2D_indices_to_1D(i,g,Ng,Nc)]*y[convert_2D_indices_to_1D(i,g,Ng,Nc)]*y[convert_2D_indices_to_1D(i,g,Ng,Nc)];}
				}
				double phi = sum1+sum2;
	
				if (phi > max_Phi) {
					max_Phi = phi;
					maxPhi_column = j;
				}
			} else {
				// Calculate Theta
				double sum = 0.0;
				int e;
				for (e = 0; e < *Ng; e++) {
					sum = sum + w[convert_2D_indices_to_1D(e,j,Ng,Nc)]*y[convert_2D_indices_to_1D(e,j,Ng,Nc)]*y[convert_2D_indices_to_1D(e,j,Ng,Nc)];
				}
				double theta = sum;
				if (theta > max_Theta) {
					max_Theta = theta;
					maxTheta_column = j;
				}
			}
		}
	}
	if (exists_zero) {
		int i;
		for (i = 0; i < *Ng; i++) {
			if (w[convert_2D_indices_to_1D(i,maxPhi_column,Ng,Nc)] > close_enough) {
				a[i] = y[convert_2D_indices_to_1D(i,maxPhi_column,Ng,Nc)];
			} else {
				// Calculate Beta
				double num = 0.0;
				double denom = 0.0;
				int e,g;
				for (e = 0; e < *Ng; e++) {
					for (g = 0; g < *Nc; g ++) {
						if (e != i && g !=  maxTheta_column) {
							num = num + w[convert_2D_indices_to_1D(e,g,Ng,Nc)]*y[convert_2D_indices_to_1D(e,maxPhi_column,Ng,Nc)]*y[convert_2D_indices_to_1D(e,maxPhi_column,Ng,Nc)]*y[convert_2D_indices_to_1D(i,g,Ng,Nc)]*y[convert_2D_indices_to_1D(i,g,Ng,Nc)];
							denom = denom + w[convert_2D_indices_to_1D(e,g,Ng,Nc)]*y[convert_2D_indices_to_1D(e,maxPhi_column,Ng,Nc)]*y[convert_2D_indices_to_1D(i,g,Ng,Nc)]*y[convert_2D_indices_to_1D(e,g,Ng,Nc)];
						}
					}
				}
				a[i] = num/denom;
			}
		}
		
	} else {
		int i;
		for (i = 0; i < *Ng; i++) {
			a[i] = y[convert_2D_indices_to_1D(i,maxTheta_column,Ng,Nc)];
		}
	}
}

int converge(double* a, double* last_a, int* Ng, double* e) {
	double distance = 0.0;
	int i;
	for (i = 0; i < *Ng; i++) {
		// euclidean distance
		distance = distance + (a[i]-last_a[i])*(a[i]-last_a[i]);
	}
	if (sqrt(distance) < *e) {
		return(1);
	} else {
		return(0);
	}
}

void calc_b (double* y, double* w, int* Ng, int* Nc, double* a, double* b ) {
	int i,j;
	for (j = 0; j < *Nc; j++) {
		double num = 0.0;
		double denom = 0.0;
		for (i = 0; i < *Ng; i++) {
			num = num + w[convert_2D_indices_to_1D(i,j,Ng,Nc)]*a[i]*y[convert_2D_indices_to_1D(i,j,Ng,Nc)];
			denom = denom + w[convert_2D_indices_to_1D(i,j,Ng,Nc)]*a[i]*a[i];
		}
		b[j] = num/denom;
	}
}

void calc_a (double* y, double* w, int* Ng, int* Nc, double* b, double* a ) {
	int i,j;
	for (i = 0; i < *Ng; i++) {
		double num = 0.0;
		double denom = 0.0;
		for (j = 0; j < *Nc; j++) {
			num = num + w[convert_2D_indices_to_1D(i,j,Ng,Nc)]*b[j]*y[convert_2D_indices_to_1D(i,j,Ng,Nc)];
			denom = denom + w[convert_2D_indices_to_1D(i,j,Ng,Nc)]*b[j]*b[j];
		}
		a[i] = num/denom;
	}
}

/* Cannot pass matrices to C from R instead matrix = long array of col1,col2,col3,etc...*/
void PCA_WGT (double* y, double* w, int* Ng, int* Nc, double* e, int* Nf, int* max_steps, double* A, double* B) {
	int f = 0;
	double* a = malloc((*Ng) * sizeof (double));
	double* last_a = malloc((*Ng) * sizeof (double));
	double* b = malloc((*Nc) * sizeof (double));
	double* last_b = malloc((*Nc) * sizeof (double));
	int used_steps = 0;
	for (f = 0; f < *Nf; f++) {
		//Calculate each component
		initialize(y,w,Ng,Nc, a);
		int count = 0;
		while (!converge(a,last_a,Ng,e) || !converge(b,last_b,Nc,e)) {
			calc_b(y,w,Ng,Nc,a, b);
			// swap pointers a
			double* tmp = last_a;
			last_a = a;
			a = tmp;
			calc_a(y,w,Ng,Nc,b, a);
			// swap pointers b
			tmp = last_b;
			last_b = b;
			b = tmp;
			count++;
			if (count > used_steps) {used_steps = count;}
			if (count > *max_steps) {
				error("Algorithm did not converge within specified number of steps.\n");
			}
		}
		// Store result
		int i,j;
		for (i = 0; i < *Ng; i++) {
			A[convert_2D_indices_to_1D(i,f,Ng,Nf)] = a[i];
		}
		for (j = 0; j < *Nc; j++) {
			B[convert_2D_indices_to_1D(f,j,Nf,Nc)] = b[j];
		}
		// subtract residuals
		for (i = 0; i < *Ng; i++) {
			for(j = 0; j < *Nc; j++) {
				y[convert_2D_indices_to_1D(i,j,Ng,Nc)] = y[convert_2D_indices_to_1D(i,j,Ng,Nc)]-a[i]*b[j];
			}
		}
	}
	free(a); free(last_a);
	free(b); free(last_b);
	*max_steps = used_steps;
}
