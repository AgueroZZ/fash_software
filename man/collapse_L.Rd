% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_BF.R
\name{collapse_L}
\alias{collapse_L}
\title{Collapse a Likelihood Matrix for Bayes Factor Computation}
\usage{
collapse_L(L, log = FALSE)
}
\arguments{
\item{L}{A numeric matrix representing the likelihoods. Rows correspond to datasets, and
columns correspond to mixture components (including the null component in the first column).}

\item{log}{A logical value. If \code{TRUE}, treats \code{L} as a log-likelihood matrix.}
}
\value{
A list containing:
\describe{
  \item{L_c}{A 2-column matrix where the first column corresponds to the null likelihood
  and the second column corresponds to the reweighted alternative likelihood.}
  \item{pi_hat_star}{A numeric vector of mixture weights estimated under the alternative hypothesis.}
}
}
\description{
This function collapses a likelihood matrix into a 2-column matrix, reweighting the likelihood
under the alternative hypothesis. This is done using the mix-SQP algorithm to estimate
the optimal mixture weights under the alternative hypothesis.
}
\examples{
# Example likelihood matrix (log-space)
set.seed(1)
L <- matrix(abs(rnorm(20)), nrow = 5, ncol = 4)
collapse_result <- fashr:::collapse_L(L, log = FALSE)
print(collapse_result$L_c)

}
\keyword{internal}
