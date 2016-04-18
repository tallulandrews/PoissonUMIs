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

PoisUMI_calc_difference <- function(expr_mat) {
	nrows = length(expr_mat[,1])
	ncols = length(expr_mat[1,])
	difference = matrix(rep(0.0, times =nrows*ncols), nrow=nrows, ncol=ncols)
	out <- .C("calc_difference",expr_mat=as.integer(expr_mat), nc = as.integer(ncols), ng = as.integer(nrows), as.double(difference));
	output = matrix(out$difference, nrow=nrows, ncol=ncols)
	colnames(output) <- colnames(expr_mat)
	rownames(output) <- rownames(expr_mat)
	return(output)
}

