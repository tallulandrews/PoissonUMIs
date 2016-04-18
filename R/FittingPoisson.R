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
	if (sum(expr_mat < 0) >0) {stop("Expression matrix contains negative values! M3Drop requires an expression matrix that is not log-transformed.")}
	p = apply(expr_mat,1,function(x){y = x[!is.na(x)]; sum(y==0)/length(y)});
	s = rowMeans(expr_mat, na.rm=T);
	return(list(s = s,p = p));
}

PoisUMI_Fit_Basic_Poisson <- function(expr_mat) {
	vals = hidden_calc_s_and_p(expr_mat)

	p = vals$p
	s = vals$s
	LL <- function(alpha, sigma) {
		p_zero = 1 - (1-exp(-alpha*s))
		R = p - p_zero;
	        R = suppressWarnings(dnorm(R, 0, sigma, log = TRUE))
	        -sum(R)
	    }
	fit = mle2(LL, start = list(alpha=1, sigma=0.1))
	return(list(s=s,p_obs=p,p_exp=exp(-fit@coef[1]*s), alpha=fit@coef[1]))
}

PoisUMI_Poisson_Account_For_Depth <- function(expr_mat, dimension="genes", scale_factor=1) {
	# Check input
	dimension = as.character(dimension)

	# Get values for later
	tis = colSums(expr_mat) # Total molecules/cell
	sj = rowSums(expr_mat)/length(expr_mat[1,]) # Mean expression each gene
	nc = length(expr_mat[1,]) # Number of cells
	ng = length(expr_mat[,1]) # Number of genes
	total = sum(tis) # Total molecules sampled

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
		return(gis_calculated, var=gis_variance);
	}
}

PoisUMI_Poisson_Account_For_Depth_DE <- function(expr_mat) {
	sj = rowSums(expr_mat)/length(expr_mat[1,]) # Mean expression each gene
	djs = rowSums(expr_mat == 0) # Observed Dropouts per cell
	nc = length(expr_mat[1,])
	drop_rate_obs = djs/nc
	drop_rate_obs_err = sqrt(drop_rate_obs*(1-drop_rate_obs)/nc)

	Expected = PoisUMI_Poisson_Account_For_Depth(expr_mat, dimension="genes")
	drop_rate_exp = Expected$vals/nc
	drop_rate_exp_err = sqrt(Expected$var/nc^2)
	difference = drop_rate_obs-drop_rate_exp
	combined_err = drop_rate_exp_err
	Zed = difference/combined_err
	pvalue = pnorm(Zed, lower.tail = FALSE)
	Expected$pvalue = pvalue
	return(list(s=sj, p_obs=djs, p_exp=Expected$vals, p_exp_var=Expected$var, pvalue=pvalue));
}


# Scale & different depth

PoisUMI_Fit_Full_Poisson <- function(expr_mat) {
	# Get values for later
	tis = colSums(expr_mat) # Total molecules/cell
	sj = rowSums(expr_mat)/length(expr_mat[1,]) # Mean expression each gene
	nc = length(expr_mat[1,]) # Number of cells
	ng = length(expr_mat[,1]) # Number of genes
	total = sum(tis) # Total molecules sampled
	gis = colSums(expr_mat==0)
	djs = rowSums(expr_mat==0)

	LL <- function(alpha) {
		R = 0
		for (j in 1:ng) {
			p = exp(-tis*sj[j]*nc*alpha/total)
			res = djs[j] - sum(p);
			sigma = sqrt(sum(0.5*(1-0.5))); # Prevent sigma being too small & causing errors.
			if (is.na(dnorm(res, 0, sigma, log=TRUE))) {print(j);}
			R = R + dnorm(res, 0, sigma, log=TRUE);
		}
		R = R/ng

		R2 = 0
		for (i in 1:nc) {
			p = exp(-tis[i]*(sj*nc*alpha/total))
			res = gis[i]-sum(p)
			sigma = sqrt(sum(p*(1-p)));
			R2 = R2+dnorm(res,0, sigma, log=TRUE)
		}
		R2 = R2/nc
	        -(R+R2)
	    }
	fit_basic = PoisUMI_Fit_Basic_Poisson(expr_mat)
	fit = mle2(LL, start = list(alpha=fit_basic$alpha))
	scale_factor = fit@coef[1]
	predictions = PoisUMI_Poisson_Account_For_Depth(expr_mat, dimension="genes", scale_factor=scale_factor)
	return(list(s=sj, p_obs = djs, p_exp=predictions$vals, p_exp_var=predictions$var, alpha=scale_factor, alpha_basic=fit_basic$alpha));
}

PoisUMI_Full_Poisson_DE <- function(expr_mat) {
	djs = rowSums(expr_mat == 0) # Observed Dropouts per cell
	nc = length(expr_mat[1,]) # Number of cells
	drop_rate_obs = djs/nc
	drop_rate_obs_err = sqrt(drop_rate_obs*(1-drop_rate_obs)/nc)

	Expected = PoisUMI_Fit_Full_Poisson(expr_mat)
	drop_rate_exp = Expected$p_exp/nc
	drop_rate_exp_err = sqrt(Expected$p_exp_var/nc^2)
	difference = drop_rate_obs-drop_rate_exp
	combined_err = drop_rate_exp_err
	Zed = difference/combined_err
	pvalue = pnorm(Zed, lower.tail = FALSE)
	Expected$pvalue = pvalue;
	return(Expected)
}
