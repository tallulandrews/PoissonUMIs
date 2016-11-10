\name{PoisUMI_Calc_Weighted_Distances}
\alias{PoisUMI_Calc_Weighted_Distances}
\title{Calculate Weighted Distances}
\description{
	Calculates a weighted distance between all pairs of cells.
}
\usage{PoisUMI_Calc_Weighted_Distances(norm_mat, weight_mat, exponent=2)}
\arguments{
  \item{norm_mat}{a numeric matrix of normalized counts, columns = samples, rows = genes, see Details.}
  \item{weight_mat}{a numeric matrix of weights of equal dimension to norm_mat.}
  \item{exponent}{exponent or norm (p) of the distance, 2 = Euclidean.}
}
\details{
	This function calculates the pair-wise weighted distance between all columns of norm_mat (cells) as
	the p-th root of the sum of absolute difference of expression for each gene take to the p-th power. 
	Where p is the exponent or norm of the distance, and exponent = 1 corresponds to Manhattan distance 
	and exponent = 2 corresponds to Euclidean distance.

	\deqn{D_{AB} = \sqrt[\frac{1}{p}]{ \sum |x_A - x_B|^p }}{D = ( sum(|x_A - x_B|^p) )^1/p}
}
\value{A matrix of pair-wise distances between samples.}
\examples{
  library(M3DExampleData)
  fit <- PoisUMI_Fit_Full_Poisson(Mmus_example_list$data)
  weights <- PoisUMI_Calc_Weighted_Distances(Mmus_example_list$data, fit$lambdas, fit$alpha)
  norm <- PoisUMI_Normalize_Data(Mmus_example_list$data, fit)
  dist <- PoisUMI_Calc_Weighted_Distances(norm, weights, 1)
}
\keyword{Poisson, distance, weighted}
