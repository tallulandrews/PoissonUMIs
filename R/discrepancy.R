PoisUMI_calc_discrepancy <- function(expr_mat) {
	nrows = length(expr_mat[,1])
	ncols = length(expr_mat[1,])
	discrepancy = matrix(rep(0.0, times =nrows*ncols), nrow=rows, ncol=cols)
	out <- .C("discrepancy",expr_mat=as.integer(expr_mat), nc = as.integer(ncols), ng = as.integer(nrows), as.double(discrepancy));
	output = matrix(out$discrepancy, nrow=nrows, ncol=ncols)
	colnames(output) <- colnames(expr_mat)
	rownames(output) <- rownames(expr_mat)
	return(output)
}

