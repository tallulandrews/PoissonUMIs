#include <math.h>
#define ITMAX 100 //max iterations
#define EPS 1.0e-3 //accuracy
#define FPMIN 1.0e-30 //Number near smallest possilbe floating-point number.

int convert_2D_indices_to_1D (int i, int j, int* nrow, int* ncol) {
	return(j* (*nrow) + i);
}

/* Returns ln(gamma(xx)) for xx > 0  with error < 2*10^-10 
From: Numerical Recipies in C.
*/
float gammln(float xx) {
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

float dPoisson (int x, float lambda) {
	if (x == 0) {
		return(exp(-lambda));
	}
	double lnp = x*log(lambda) - lambda - (gammln(x+1.0));
	return(exp(lnp));
}

/* Poisson cdf from Numerical Recipes in C */
float pPoisson (int k, float lambda) {
	void gcf(float* gammcf, float a, float x, float* gln);
	void gsr(float* gamser, float a, float x, float* gln);
	void nrerror(char error_text[]);
	float gamser, gammcf, gln;

	if (lambda < 0.0 || k < 0) nerror("Invalid arguments in routine gammq");
	if (k == 0) {
		return(dPoisson(k, lambda));
	}
	if (k < 25) {
		float sum = 0.0;
		int i;
		for (i = 0; i <=k; i++) {
			sum+=dPoisson(i, lambda);
		}
		return(sum);
	}
	if (lambda < (k + 1.0)) {
		gser(&gamser, (float) k, lambda, &gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf, (float) k, lambda, &gln);
		return gammcf;
	}
}

void gser(float* gamser, float a, float x, float* gln) {
	float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
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

void gcf(float* gammcf, float a, float x, float* gln) {
	float gammln(float xx);
	void nrerror(char error_text[]);
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FMIN) d=FPMIN;
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

