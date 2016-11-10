\name{PoisUMI_Normalize_Data}
\alias{PoisUMI_Normalize_Data}
\title{Normalize Data}
\description{
	Calculates weights based on information contained in each expression value (i.e. difference from null expectation).
}
\usage{PoisUMI_Normalize_Data(expr_mat, fit=NA)}
\arguments{
  \item{expr_mat}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{fit}{output from PoisUMI_Fit_Full_Poisson.}
}
\details{The Poisson model considers each observation to be drawn from a Poisson distribution with parameter:
	 \deqn{\lambda_{ij} = m_i*m_j*T*\alpha}{lambda_ij = m_i * m_j * T * alpha}
	where \deqn{m_i}{m_i} and \deqn{m_j}{m_j} are the proportion of all molecules detected that are of gene i or in cell j respectively, 
	and T is the total molecules detected across all genes and cells. \deqn{\alpha}{alpha} approximates the extent of over-counting of molecules.

	This function normalizes expression data as the number of standard deviations from the mean each observation is using the observation-specfic 
	Poisson distribution. Resulting data will contain positive, higher than expected, and negative, lower than expected, values.

	If "fit" is not provided it will be automatically calculated using PoisUMI_Fit_Full_Poisson.

}
\value{A matrix of normalized expression values.}
\examples{
  library(M3DExampleData)
  fit <- PoisUMI_Fit_Full_Poisson(Mmus_example_list$data)
  norm <- PoisUMI_Normalize_Data(Mmus_example_list$data, fit)
} 
\keyword{Poisson, normalization}

