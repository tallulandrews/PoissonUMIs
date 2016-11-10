
# Do I want to include this?
PoisUMI_Poisson_Account_For_Depth_DE <- function(expr_mat) {
	vals = hidden_calc_s_and_p(expr_mat)
	sj = vals$s # Mean expression each gene
	djs = vals$djs # Observed Dropouts per gene
	nc = vals$nc
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


PoisUMI_Full_Poisson_DE <- function(expr_mat, fit=NA) {
	# Need to incorporate error on S
	vals = hidden_calc_s_and_p(expr_mat)
	djs = vals$djs # Observed Dropouts per gene
	nc = vals$nc # Number of cells
	drop_rate_obs = djs/nc
	drop_rate_obs_err = sqrt(drop_rate_obs*(1-drop_rate_obs)/nc)

	if (!("p_exp" %in% names(fit))) {
		Expected = PoisUMI_Fit_Full_Poisson(expr_mat)
	} else {
		Expected=fit
	}
	drop_rate_exp = Expected$p_exp/nc
	drop_rate_exp_err = sqrt(Expected$p_exp_var/nc^2)
	difference = drop_rate_obs-drop_rate_exp
	combined_err = drop_rate_exp_err
	Zed = difference/combined_err
	pvalue = pnorm(Zed, lower.tail = FALSE)
	Expected$pvalue = pvalue;
	return(Expected)
}

PoisUMI_Important_Genes <- function(expr_mat, fit=NULL, method="weight") {
	if (!(method %in% c("weight","outlier"))) {
		stop(paste("Error: unknown method:", method, "\n"));
	}
	if (is.null(fit)) {
		fit = PoisUMI_Fit_Full_Poisson(expr_mat)
	}
	if (method == "weight") {
		weights = PoisUMI_Calc_Weights(expr_mat, fit$lambdas, fit$alpha)
		score = rowSums(abs(weights));
	}
	if (method == "outlier") {
		expect = PoisUMI_Full_Poisson_DE(expr_mat, fit=fit)
		score = qnorm(expect$pvalue, lower.tail=FALSE)
	}
	names(score) = rownames(expr_mat);
	return(rev(sort(score)));
}

PoisUMI_Segment_DE_Genes <- function(DE_table, fdr=0.05, tolerance=0.2) {
	ncol = length(DE_table);
	sig <- DE_table[DE_table[,ncol] < fdr,]; #Significant DE genes
	sig = sig[,-c(ncol, ncol-1)]; # remove pval & qval from table
	ncol = length(sig[1,]);
	
	for (gene in 1:length(sig[,1])) {
		vals = sig[gene,];
		ranks = rank(vals);
		max_gap = vals[which(ranks == max(ranks))] - vals[which(ranks == min(ranks))]
		for(i in 2:max(ranks)) {

			next_highest = max(ranks[ranks < i])
			gap = vals[which(ranks==i)] - vals[which(ranks == next_highest)];
			if (sum(gap < max_gap*tolerance)) {
				ranks[ranks==i] = next_highest;
			}
		}
		
		

	}
}
