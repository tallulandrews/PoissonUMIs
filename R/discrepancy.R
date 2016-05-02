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
