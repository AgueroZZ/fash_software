% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_fash.R
\name{fash}
\alias{fash}
\title{Perform Full FASH Analysis}
\usage{
fash(
  Y = NULL,
  smooth_var = NULL,
  offset = 0,
  S = NULL,
  Omega = NULL,
  data_list = NULL,
  grid = seq(0, 2, length.out = 10),
  likelihood = "gaussian",
  num_basis = 30,
  betaprec = 1e-06,
  order = 2,
  pred_step = 1,
  penalty = 1,
  num_cores = 1,
  verbose = FALSE
)
}
\arguments{
\item{Y}{Either a numeric matrix of response variables or a character string specifying the column name in `data_list` for response variables.}

\item{smooth_var}{A numeric matrix, vector, or a character string specifying the column name in `data_list` for smoothing variables.}

\item{offset}{A numeric matrix, vector, scalar, or a character string specifying the column name in `data_list` for offset variables. Defaults to 0.}

\item{S}{A numeric matrix, vector, scalar, or list representing the standard errors of `Y`. Or a character string specifying the column name in `data_list` for SD. Default is `NULL`.}

\item{Omega}{Either a list of precision matrices (one for each dataset) or a single precision matrix (shared across all datasets). Default is `NULL`.}

\item{data_list}{A list of data frames, where each data frame corresponds to a single dataset. Default is `NULL`.}

\item{grid}{A numeric vector representing the grid of PSD (Predictive Standard Deviation) values. Default is `seq(0, 2, length.out = 10)`.}

\item{likelihood}{A character string specifying the likelihood function to use. Options are `"gaussian"` and `"poisson"`. Default is `"gaussian"`.}

\item{num_basis}{An integer specifying the number of O-Spline basis functions. Default is 30.}

\item{betaprec}{A numeric value representing the precision of the fixed effects coefficients (`beta`). Default is `1e-6`.}

\item{order}{An integer specifying the order of the Integrated Wiener Process (IWP) prior. Default is 2.}

\item{pred_step}{A numeric value specifying the prediction step size. Default is 1.}

\item{penalty}{A numeric value representing the lambda value for the Dirichlet prior. Default is 1 (no additional penalty).}

\item{num_cores}{An integer specifying the number of cores to use for parallel processing. Default is 1.}

\item{verbose}{A logical value. If `TRUE`, shows progress messages and timing for each step. Default is `FALSE`.}
}
\value{
A `fash` object containing:
  \describe{
    \item{\code{prior_weights}}{Estimated prior weights for PSD values.}
    \item{\code{posterior_weights}}{Posterior weight matrix of PSD values.}
    \item{\code{psd_grid}}{PSD grid values.}
    \item{\code{lfdr}}{Local false discovery rate for each dataset.}
    \item{\code{settings}}{A list of settings used in the FASH pipeline.}
    \item{\code{fash_data}}{A structured list of data components.}
    \item{\code{L_matrix}}{Likelihood matrix used in the FASH pipeline.}
    \item{\code{eb_result}}{Empirical Bayes estimation results.}
  }
}
\description{
This function performs the full FASH pipeline, including data setup, likelihood computation,
empirical Bayes estimation, and outputs a structured `fash` object.
}
\examples{
data_list <- list(
  data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
  data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
)
grid <- seq(0, 2, length.out = 10)
result <- fash(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", grid = grid, likelihood = "poisson", verbose = TRUE)
print(result)

}
