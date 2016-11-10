#Copyright (c) 2016 Genome Research Ltd .
#Author : Tallulah Andrews <tallulandrews@gmail.com>
#This file is part of PoissonUMIs.

#PoissonUMIs is free software : you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation; either version 2 of the License, or (at your option) any later
#version.

#This program is distributed in the hope that it will be useful, but WITHOUT
#ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with
#this program . If not , see <http://www.gnu.org/licenses/>.


PoisUMI_Down_Sample_Matrix <- function(expr_mat) {
	min_lib_size = min(colSums(expr_mat))
	down_sample <- function(x) {
		prob = min_lib_size/sum(x)
		return(unlist(lapply(x, function(y) {rbinom(1,y,prob)})))
	}
	down_sampled_mat = apply(expr_mat,2,down_sample)
	return(down_sampled_mat)
}

hidden_calc_s_and_p <- function(expr_mat) {
	if (sum(expr_mat < 0) >0) {stop("Expression matrix contains negative values! PoissonUMIs requires an expression matrix that is not log-transformed.")}
	p = apply(expr_mat,1,function(x){y = x[!is.na(x)]; sum(y==0)/length(y)});
	s = rowMeans(expr_mat, na.rm=T);

	tis = colSums(expr_mat, na.rm=T) # Total molecules/cell
	djs = rowSums(expr_mat == 0, na.rm=T) # Observed Dropouts per gene
	gis = colSums(expr_mat == 0, na.rm=T) # Observed Dropouts per cell
	nc = length(expr_mat[1,]) # Number of cells
	ng = length(expr_mat[,1]) # Number of genes
	total = sum(tis, na.rm=T) # Total molecules sampled
	return(list(s = s,p = p, tis=tis, total=total,nc=nc,ng=ng, djs=djs, gis=gis));
}

PoisUMI_Fit_Basic_Poisson <- function(expr_mat, sigma_init=NA) {
	vals = hidden_calc_s_and_p(expr_mat)
	nc = length(expr_mat[1,]) # Number of cells
        ng = length(expr_mat[,1]) # Number of genes


	p = vals$p
	s = vals$s
	if (is.na(sigma_init)) {sigma_init = mean(p*(1-p))}
	LL <- function(alpha, sigma) {
		if (alpha <= 0) { return(nc*ng*1000000000)} #alpha should always be positive
		if (sigma <= 0) { return(nc*ng*1000000000)} #sigma should always be positive
		p_zero = 1 - (1-exp(-alpha*s))
		if (sum(is.na(p_zero) | p_zero < 0 | p_zero > 1) > 0) {
			stop("p is unacceptable value (1)")
		}
		R = p - p_zero;
	        R = dnorm(R, 0, sigma, log = TRUE)
	        -sum(R)
	    }
	fit = mle2(LL, start = list(alpha=1, sigma=sigma_init))
	return(list(s=s,p_obs=p,p_exp=exp(-fit@coef[1]*s), alpha=fit@coef[1]))
}

PoisUMI_Poisson_Account_For_Depth <- function(expr_mat, dimension="genes", scale_factor=1) {
	# Check input
	dimension = as.character(dimension)
	vals = hidden_calc_s_and_p(expr_mat)

	# Get values for later
	tis = vals$tis # Total molecules/cell
	sj = vals$s # Mean expression each gene
	nc = vals$nc # Number of cells
	ng = vals$ng # Number of genes
	total = vals$total # Total molecules sampled

	if (dimension == "1" | dimension == "genes" | dimension == "rows") {
		djs_calculated = vector(length=ng)
		djs_variance = vector(length=ng)
		for (j in 1:ng) {
			p = exp(-tis*sj[j]*scale_factor*nc/total)
			djs_calculated[j] = sum(p)
			djs_variance[j] = sum(p*(1-p))
		}
		return(list(vals=djs_calculated, var=djs_variance))
	} else if (dimension == "2" | dimension == "cells" | dimension == "cols") {
		gis_calculated = vector(length=nc)
		gis_variance = vector(length=nc)
		for (i in 1:nc) {
			p = exp(-tis[i]*(sj*scale_factor*nc/total))
			gis_calculated[i] = sum(p)
			gis_variance[i] = sum(p*(1-p))
		}
		return(list(vals=gis_calculated, var=gis_variance));
	}
}



# Scale & different depth

PoisUMI_Fit_Full_Poisson <- function(expr_mat) {
	vals = hidden_calc_s_and_p(expr_mat)
	# Get values for later
	tis = vals$tis # Total molecules/cell
	sj = vals$s # Mean expression each gene
	nc = vals$nc # Number of cells
	ng = vals$ng # Number of genes
	total = vals$total # Total molecules sampled
	gis = vals$gis # Number dropouts per cell
	djs = vals$djs # Number dropouts per gene

	LL <- function(alpha) {
		if (alpha <= 0) { return(nc*ng*1000000000)} #alpha should always be positive
		R = 0;
		denom = ng;
		for (j in 1:ng) {
			p = exp(-tis*sj[j]*nc*alpha/total)
			if (sum(is.na(p) | p < 0 | p > 1) > 0) {
				stop("p is unacceptable value (1)")
			}
			res = djs[j] - sum(p);
			sigma = sqrt(sum(p*(1-p))); 
			sigma_obs =  sqrt(nc*djs[j]/nc*(1-djs[j]/nc))
			sigma = max(sigma, sigma_obs)
			if (sigma == 0 & res == 0){R = R+1; next;} 
			if (sigma == 0 & abs(res) > 0){denom = denom-1; next;} 
#			if (is.na(dnorm(res, 0, sigma, log=TRUE))) {print(j);}
#			if (!is.finite(dnorm(res, 0, sigma, log=TRUE))) {print(j);}
			R = R + dnorm(res, 0, sigma, log=TRUE);
		}
		R = R/denom

		R2 = 0
		denom2 = nc
		for (i in 1:nc) {
			p = exp(-tis[i]*(sj*nc*alpha/total))
			if (sum(is.na(p) | p < 0 | p > 1) > 0) {
				stop("p is unacceptable value (2)")
			}
			res = gis[i]-sum(p)
			sigma = sqrt(sum(p*(1-p)));
			sigma_obs = sqrt(ng*gis[i]/ng*(1-gis[i]/ng))
			sigma = max(sigma, sigma_obs)
			if (sigma == 0 & res == 0){R2=R2+1; next;} 
			if (sigma == 0 & abs(res) > 0){denom2 = denom2-1; next;} 
			R2 = R2+dnorm(res,0, sigma, log=TRUE)
		}
		R2 = R2/denom2
	        -(R+R2)
	    }
	# Fit poisson without accounting for cell-specific depth to find approx scaling
	fit_basic = PoisUMI_Fit_Basic_Poisson(expr_mat)
	# Add robustness against failure to fit
	res = fit_basic$p_obs - fit_basic$p_exp;
	if (sum(res < 0)/length(res) > 0.75 | sum(res > 0)/length(res) > 0.75) {
		print("First fitting failed attempting again")
		fit_basic = PoisUMI_Fit_Basic_Poisson(expr_mat, sigma_init=0.01)
	}
	res = fit_basic$p_obs - fit_basic$p_exp;
	if (sum(res < 0)/length(res) > 0.75 | sum(res > 0)/length(res) > 0.75) {
		print("Second fitting failed attempting again")
		fit_basic = PoisUMI_Fit_Basic_Poisson(expr_mat, sigma_init=0.5)
	}

	# Fit the full model starting at fit above b/c this fitting is slow
	fit = mle2(LL, start = list(alpha=fit_basic$alpha))
	scale_factor = fit@coef[1]
	predictions = PoisUMI_Poisson_Account_For_Depth(expr_mat, dimension="genes", scale_factor=scale_factor)
	
	# Calculate lambda for each observation
	lambdas = (sj %*% t(tis))*(nc*scale_factor/total);

	# Calculate log-normal error of alpha
	gene_specific_alpha <- function(gene_index) {
		LLgs <- function(alpha) {
			if (alpha <= 0) { return(nc*ng*1000000000)} #alpha should always be positive
			p = exp(-tis*sj[gene_index]*nc*alpha/total)
			if (sum(is.na(p) | p < 0 | p > 1) > 0) {
				stop("p is unacceptable value (1)")
			}
			sigma = sqrt(sum(p*(1-p)));
			sigma_obs =  sqrt(nc*djs[gene_index]/nc*(1-djs[gene_index]/nc))
			sigma = max(sigma, sigma_obs)
			res = djs[gene_index] - sum(p);
			if (sigma == 0 & res == 0){return(0)} 
			if (sigma == 0 & abs(res) > 0){return(10000000000)} 
			-1*dnorm(res, 0, sigma, log=TRUE);
		}
		fit = mle2(LLgs, start = list(alpha=scale_factor))
		return(fit@coef[1]);
	}
	alphas = sapply(1:length(sj), function(x){suppressWarnings(gene_specific_alpha(x))});
	alpha_err = sqrt(var(log(alphas))/length(sj)) # Is this correct? - Methinks not

	return(list(s=sj, p_obs = djs, p_exp=predictions$vals, p_exp_var=predictions$var, alpha=scale_factor, alpha_basic=fit_basic$alpha, alpha_err=alpha_err, lambdas = lambdas));
}

PoisUMI_calc_lambdas <- function(expr_mat) {
	fit = PoisUMI_Fit_Full_Poisson(expr_mat);
	return(fit$lambdas);
}

