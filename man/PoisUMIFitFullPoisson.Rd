\name{PoisUMI_Fit_Full_Poisson}
\alias{PoisUMI_Fit_Full_Poisson}
\title{Fits the Poisson UMI model}
\description{
	Fits the Poisson model's scaling factor (alpha) using the observed relationship between mean-expression and dropout rate.
}
\usage{PoisUMI_Fit_Full_Poisson(expr_mat)}
\arguments{
  \item{expr_mat}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
}
\details{The Poisson model considers each observation to be drawn from a Poisson distribution with parameter:
	 \deqn{\lambda_{ij} = m_i*m_j*T*\alpha}{lambda_ij = m_i * m_j * T * alpha}
	where \deqn{m_i}{m_i} and \deqn{m_j}{m_j} are the proportion of all molecules detected that are of gene i or in cell j respectively, 
	and T is the total molecules detected across all genes and cells. 

	This function fits the scaling parameter \deqn{\alpha}{alpha} which approximates the extent of over-counting of molecules using the observed
	relationship between mean-expression and dropout rate. Fitting is performed in two steps using maximum-likelihood estimation implemented by 
	mle2 from the statmod package. First, a simplified model which excludes the \deqn{m_j}{m_j} factor is fit using multiple starting points in 
	case of poor convergence. Second, the estimated \deqn{\alpha}{alpha} is used as a starting point to fit the full model. 

	This procedure improves stability and speed of fitting.
}
\value{List of output containing:
	s = observed mean expression of each gene (m_i)
	p_obs = observed gene-specific dropout rate
	p_exp = model-based gene-specific dropout rate
	p_exp_var = model-based variance of gene-specific dropout rate
	alpha = final fit scaling parameter
	alpha_basic = initial fit scaling parameter from simplified model
	lambdas = matrix of calculated lambda for each expression value (same dimensions as input matrix).
}
\examples{
  library(M3DExampleData)
  fit <- PoisUMI_Fit_Full_Poisson(Mmus_example_list$data)
} 
\keyword{Poisson, model fitting, maximum likelihood}
\seealso{mle2}
