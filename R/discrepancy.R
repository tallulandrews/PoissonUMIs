PoisUMI_calc_discrepancy <- function(expr_mat) {
	nrows = length(expr_mat[,1])
	ncols = length(expr_mat[1,])
	discrepancy = matrix(rep(0.0, times =nrows*ncols), nrow=nrows, ncol=ncols)
	out <- .C("calc_discrepancy",expr_mat=as.integer(expr_mat), nc = as.integer(ncols), ng = as.integer(nrows), as.double(discrepancy));
	output = matrix(out$discrepancy, nrow=nrows, ncol=ncols)
	colnames(output) <- colnames(expr_mat)
	rownames(output) <- rownames(expr_mat)
	return(output)
}

PoisUMI_pearson_cor <- function(expr_mat, lambdas=NA) {
	nrows = length(expr_mat[,1])
	ncols = length(expr_mat[1,])
	if (is.na(lambdas)) {
		lambdas = PoisUMI_calc_lambdas(expr_mat);
	}
	cor_mat = matrix(rep(0.0, times = ncols*ncols), nrow=ncols, ncol=ncols)
	out <- .C("calc_difference",expr_mat=as.integer(expr_mat), nc = as.integer(ncols), ng = as.integer(nrows), as.double(lambdas), as.double(cor_mat));
	output = matrix(out$cor_mat, nrow=nrows, ncol=ncols)
	colnames(output) <- colnames(expr_mat)
	rownames(output) <- colnames(expr_mat)
	return(output)
}

#data = read.table("/lustre/scratch108/compgen/team218/TA/Bergiers_Wafergen/Bergiers_Wafergen_Transcriptome_DirAdj_wMulti_FeatureCounts.txt", header=T)
#fit = PoisUMI_Fit_Full_Poisson(as.matrix(data))
#library(skellam)
#PoisUMI_distance(data[,1],data[,600],lambdas[,1],lambdas[,600])
#PoisUMI_distance(data_keep[2,],data_keep[1,],lambda_keep[2,],lambda_keep[1,])


# This is symetrical, D(1,2)=D(2,1), and in range [0,Inf], and D(1,1) == 0 
# However this is slow for practical applications. 
PoisUMI_distance <- function(Obs1, Obs2, Lambda1, Lambda2, alpha) {
	Obs1 = alpha*Obs1
	Obs2 = alpha*Obs2
	# Calculate p-value for each pair of elements
	difference = unlist(Obs1-Obs2);
	thing = cbind(abs(difference), Lambda1, Lambda2)
	density_of_tails <- function(val, l1, l2) {
		if (val == 0) {
			return(1)
		} else {
			p_upper = pskellam(floor(val),lambda1=l1,lambda2=l2, lower.tail=F)
			p_lower = pskellam(-ceiling(val),lambda1=l1,lambda2=l2, lower.tail=T)
			return(p_upper+p_lower)
		}
	}
	ps = apply(thing,1, function(x){density_of_tails(x[1],l1=x[2],l2=x[3])})
	# Combine p-values using Fisher's method
	chi_sq = -2*sum(log(ps))
	df = 2*length(ps)
	pchisq(chi_sq, df=df, lower.tail=FALSE)
	return(chi_sq/df)
}
