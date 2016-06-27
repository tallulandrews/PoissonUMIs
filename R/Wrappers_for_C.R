PoisUMI_Calc_Weights <- function(expr_mat, lambdas, alpha) {
	# Check Input
	if (!is.matrix(expr_mat)) {
		expr_mat = as.matrix(expr_mat);
	}
	if (!is.matrix(lambdas)) {
		lambdas = as.matrix(lambdas);
	}

	if (!identical(dim(expr_mat), dim(lambdas))) {
		stop("Error: Expression matrix and lambdas are not of the same dimension");
	}
	Nrow = dim(expr_mat)[1];
	Ncol = dim(expr_mat)[2];
	initialized_output = matrix(rep(0, times=Nrow*Ncol), nrow=Nrow, ncol=Ncol);
	# Run command
	out <- .C( "calc_weights", val_mat=as.double(expr_mat*alpha), lambda_mat=as.double(lambdas), weight_mat=as.double(intialized_output), ncol=as.integer(Ncol), nrow=as.integer(Nrow) );

	# Clean up output
	final_output = matrix(out$weight_mat, nrow=Nrow, ncol=Ncol);
	colnames(final_output) = colnames(expr_mat);
	rownames(final_output) = rownames(expr_mat);
	return(final_output);
}

PoisUMI_Calc_Weighted_Distances  <- function(norm_mat, weight_mat, exponent=2) {
	# Check Input
	if (!is.matrix(norm_mat)) {
		norm_mat = as.matrix(norm_mat);
	}
	if (!is.matrix(weight_mat)) {
		weight_mat = as.matrix(weight_mat);
	}

	if (!identical(dim(norm_mat), dim(weight_mat))) {
		stop("Error: Expression matrix and lambdas are not of the same dimension");
	}
	Nrow = dim(norm_mat)[1];
	Ncol = dim(norm_mat)[2];
	initialized_output = matrix(rep(0, times=Nrow*Nrow), nrow=Nrow, ncol=Nrow);
	# Run command
	out <- .C( "distance_wgt", y=as.double(norm_mat), w=as.double(weight_mat), nrow=as.integer(Nrow), ncol=as.integer(Ncol), exponent = as.double(exponent), out = as.double(initialized_output) );

	# Clean up output
	final_output = matrix(out$out, nrow=Nrow, ncol=Nrow);
	colnames(final_output) = rownames(norm_mat);
	rownames(final_output) = rownames(norm_mat);
	return(final_output);
}
