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
	initialized_output <- matrix(rep(0, times=Nrow*Ncol), nrow=Nrow, ncol=Ncol);
	# Run command
	out <- .C( "calc_weights", val_mat=as.double(expr_mat*alpha), lambda_mat=as.double(lambdas), weight_mat=as.double(initialized_output), ncol=as.integer(Ncol), nrow=as.integer(Nrow) );

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

	norm_mat <- t(norm_mat)
	weight_mat <- t(weight_mat)

	if (!identical(dim(norm_mat), dim(weight_mat))) {
		stop("Error: Normalized expression matrix and weights are not of the same dimension");
	}
	Nrow = dim(norm_mat)[1];
	Ncol = dim(norm_mat)[2];
	initialized_output <- matrix(rep(0, times=Nrow*Nrow), nrow=Nrow, ncol=Nrow);
	# Run command
	out <- .C( "distance_wgt", y=as.double(norm_mat), w=as.double(weight_mat), nrow=as.integer(Nrow), ncol=as.integer(Ncol), exponent = as.double(exponent), out = as.double(initialized_output) );

	# Clean up output
	final_output = matrix(out$out, nrow=Nrow, ncol=Nrow);
	colnames(final_output) = rownames(norm_mat);
	rownames(final_output) = rownames(norm_mat);
	return(final_output);
}


PoisUMI_Calc_Weighted_Correlations  <- function(norm_mat, weight_mat) {
	# Check Input
	if (!is.matrix(norm_mat)) {
		norm_mat = as.matrix(norm_mat);
	}
	if (!is.matrix(weight_mat)) {
		weight_mat = as.matrix(weight_mat);
	}

	if (!identical(dim(norm_mat), dim(weight_mat))) {
		stop("Error: Normalized expression matrix and weights are not of the same dimension");
	}
	Nrow = dim(norm_mat)[1];
	Ncol = dim(norm_mat)[2];
	initialized_output <- matrix(rep(0, times=Nrow*Nrow), nrow=Nrow, ncol=Nrow);
	# Run command
	out <- .C( "cor_matrix_weighted", y=as.double(norm_mat), w=as.double(weight_mat), nrow=as.integer(Nrow), ncol=as.integer(Ncol), out = as.double(initialized_output) );

	# Clean up output
	final_output = matrix(out$out, nrow=Nrow, ncol=Nrow);
	colnames(final_output) = rownames(norm_mat);
	rownames(final_output) = rownames(norm_mat);
	return(final_output);
}

PoisUMI_Calc_Weighted_PCA  <- function(norm_mat, weight_mat, Nf=2, e=0.001, max_steps=10000) {
	if (!is.matrix(norm_mat)) {
		norm_mat = as.matrix(norm_mat);
	}
	if (!is.matrix(weight_mat)) {
		weight_mat = as.matrix(weight_mat);
	}

	if (!identical(dim(norm_mat), dim(weight_mat))) {
		stop("Error: Normalized expression matrix and weights are not of the same dimension");
	}
	Nrow = dim(norm_mat)[1];
	Ncol = dim(norm_mat)[2];

	A = matrix(rep(0, times = Nrow*Nf), nrow = Nrow, ncol = Nf) #Scores
	B = matrix(rep(0, times = Ncol*Nf), nrow = Nf, ncol = Ncol) #Loadings
	out <- .C( "PCA_WGT", y=as.double(norm_mat), w=as.double(weight_mat), Ng=as.integer(Nrow), Nc=as.integer(Ncol), e=as.double(e), Nf = as.integer(Nf), max_steps=as.integer(max_steps), A=as.double(A), B=as.double(B) )

	# Format output
	out$y = matrix(out$y, nrow=Nrow, ncol=Ncol)
	out$w = matrix(out$w, nrow=Nrow, ncol=Ncol)
	out$A = matrix(out$A, nrow=Nrow)
	out$B = matrix(out$B, ncol=Ncol)
	colnames(out$y) = colnames(norm_mat)
	rownames(out$y) = rownames(norm_mat)
	colnames(out$B) = colnames(norm_mat)
	rownames(out$B) = paste("PC", 1:length(out$A[1,]), sep="")
	colnames(out$A) = paste("PC", 1:length(out$B[,1]), sep="")
	rownames(out$A) = rownames(norm_mat)
	out$used_steps = out$max_steps
	out$max_steps = max_steps
	names(out)[which(names(out) == "A")] = "scores";
	names(out)[which(names(out) == "B")] = "loadings";

	return(out)
}


PoisUMI_Differential_Expression <- function(expr_mat, lambdas, alpha, labels) {
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
	if (length(labels) != Ncol) {
		stop("Error: must provide exactly one label per column of expression matrix"); 
	}
	numeric_labels = as.numeric(factor(labels))-1;
	nlabels = length(unique(labels));

	initialized_output <- matrix(rep(0, times=Nrow*(nlabels+1)), nrow=Nrow, ncol=(nlabels+1));
	# Run command
	out <- .C( "diff_expression", y=as.double(expr_mat*alpha), w=as.double(lambdas), nrow=as.integer(Nrow), ncol=as.integer(Ncol), labels=as.integer(numeric_labels), nlabels = as.integer(nlabels), out = as.double(initialized_output) );

	# Clean up output
	final_output = matrix(out$out, nrow=Nrow, ncol=(nlabels+1));
	colnames(final_output) = c(levels(factor(labels)),"pvalue");
	rownames(final_output) = rownames(norm_mat);
	return(final_output);
}
