\name{PoisUMI_Calc_Weights}
\alias{PoisUMI_Calc_Weights}
\title{Fits the Poisson UMI model}
\description{
	Calculates weights based on information contained in each expression value (i.e. difference from null expectation).
}
\usage{PoisUMI_Calc_Weights(expr_mat)}
\arguments{
  \item{expr_mat}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{lambdas}{a numeric matrix of estimated lambdas of equal dimension to expr_mat, see Details.}
  \item{alpha}{scaling parameter, see Details.}
}
\details{The Poisson model considers each observation to be drawn from a Poisson distribution with parameter:
	 \deqn{\lambda_{ij} = m_i*m_j*T*\alpha}{lambda_ij = m_i * m_j * T * alpha}
	where \deqn{m_i}{m_i} and \deqn{m_j}{m_j} are the proportion of all molecules detected that are of gene i or in cell j respectively, 
	and T is the total molecules detected across all genes and cells. \deqn{\alpha}{alpha} approximates the extent of over-counting of molecules.

	This function calculates weights for each observation as the probability of observing a greater or equal deviation from the mean of the Poisson distribution.

	lambdas and alpha should be calculated using the PoisUMI_Fit_Full_Poisson function.
}
\value{A matrix of calculated weights for each expression value (same dimensions as input matrix).}
\examples{
  library(M3DExampleData)
  fit <- PoisUMI_Fit_Full_Poisson(Mmus_example_list$data)
  weights <- PoisUMI_Calc_Weights(Mmus_example_list$data, fit$lambdas, fit$alpha)
} 
\keyword{Poisson, model fitting, maximum likelihood}
\seealso{mle2}
