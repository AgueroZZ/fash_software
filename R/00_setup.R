#' Roxygen commands
#'
#' @useDynLib Gaussian_dep
#' @useDynLib Gaussian_ind
#' @useDynLib Gaussian_dep_fixed
#' @useDynLib Gaussian_ind_fixed
#' @useDynLib Poisson_ind
#' @useDynLib Poisson_ind_fixed
#'
dummy <- function() {
  return(NULL)
}


#' Setup Quantities for FASH
#'
#' This function prepares the necessary data structures for implementing the Function Adaptive SHrinkage (FASH) approach. It processes the input data, smoothing variables, offsets, standard errors, and optional precision matrices to generate a list of data frames and supporting matrices.
#'
#' @param Y Either a numeric matrix of response variables or a character string specifying the column name in `data_list` for response variables. Required if `data_list` is provided.
#' @param smooth_var A numeric matrix, vector, or a character string specifying the column name in `data_list` for smoothing variables. Required if `data_list` is provided.
#' @param offset A numeric matrix, vector, scalar, or a character string specifying the column name in `data_list` for offset variables. Defaults to 0.
#' @param S A numeric matrix, vector, scalar, or list representing the standard errors of `Y`. Or a character string specifying the column name in `data_list` for SD. If provided as a matrix, it must have the same dimensions as `Y`. If a vector or scalar, it will be expanded to match the dimensions of `Y`. Default is `NULL`.
#' @param Omega Either a list of precision matrices (one for each dataset) or a single precision matrix (shared across all datasets). Default is `NULL`.
#' @param data_list A list of data frames, where each data frame corresponds to a single dataset. Default is `NULL`.
#'
#' @return A list with the following components:
#' \item{data_list}{A list of data frames, where each data frame corresponds to a row of `Y`. Each data frame contains:
#'   \describe{
#'     \item{\code{y}}{The response variables for the corresponding row of `Y`.}
#'     \item{\code{x}}{The smoothing variables for the corresponding row of `smooth_var`.}
#'     \item{\code{offset}}{The offset values for the corresponding row of `offset`.}
#'   }
#' }
#' \item{S}{A list of standard errors, where each element corresponds to the standard errors for a single dataset.}
#' \item{Omega}{A list of precision matrices, where each element corresponds to a single dataset's precision matrix.}
#'
#' @examples
#' # Example usage with matrix input
#' Y <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' S <- c(0.5, 0.8, 1.2, 1.0, 0.9)
#' Omega <- diag(c(1, 2, 3, 4, 5))
#' result <- fash_set_data(NULL, Y, smooth_var, offset, S, Omega)
#'
#' # Example usage with data_list input
#' data_list <- list(
#'   data.frame(y = rnorm(5), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(5), x = 1:5, offset = 0)
#' )
#' S_list <- list(rep(0.5, 5), rep(0.8, 5))
#' Omega_list <- list(diag(5), diag(5))
#' result <- fash_set_data(data_list, "y", "x", "offset", S_list, Omega_list)
#'
#' @export
fash_set_data <- function(Y, smooth_var, offset = 0, S = NULL, Omega = NULL, data_list = NULL) {
  # If data_list is provided, extract variables
  if (!is.null(data_list)) {
    if (!is.list(data_list) || !all(sapply(data_list, is.data.frame))) {
      stop("data_list must be a list of data frames.")
    }
    n_datasets <- length(data_list)

    # Extract variables from data_list
    processed_data_list <- lapply(data_list, function(df) {
      # Extract response variable (Y)
      if (!Y %in% names(df)) stop(sprintf("Variable '%s' not found in data frame.", Y))
      y <- df[[Y]]

      # Extract smoothing variable (smooth_var)
      if (!smooth_var %in% names(df)) stop(sprintf("Variable '%s' not found in data frame.", smooth_var))
      x <- df[[smooth_var]]

      # Extract or process offset
      if (is.character(offset) && offset %in% names(df)) {
        offset_value <- df[[offset]]
      } else {
        # Handle offset as matrix/vector/scalar
        n <- length(y)
        if (is.matrix(offset)) {
          if (nrow(offset) != n_datasets || ncol(offset) != ncol(Y)) {
            stop("offset must have the same dimensions as Y if it is a matrix.")
          }
        } else if (length(offset) == 1) {
          offset <- rep(offset, n)
        } else {
          stop("When using data_list, offset must be a valid column name or a scalar value.")
        }
        offset_value <- offset
      }

      # Return processed data frame
      data.frame(y = y, x = x, offset = offset_value)
    })

  } else {
    # If data_list is NULL, construct it from Y, smooth_var, and offset
    n_datasets <- nrow(Y)

    # Validate and process smooth_var
    if (is.matrix(smooth_var)) {
      if (nrow(smooth_var) != n_datasets || ncol(smooth_var) != ncol(Y)) {
        stop("smooth_var must have the same dimensions as Y if it is a matrix.")
      }
    } else if (is.vector(smooth_var) && length(smooth_var) == ncol(Y)) {
      smooth_var <- matrix(rep(smooth_var, n_datasets), nrow = n_datasets, byrow = TRUE)
    } else {
      stop("smooth_var must be a matrix or a vector matching the number of columns in Y.")
    }

    # Validate and process offset
    if (is.matrix(offset)) {
      if (nrow(offset) != n_datasets || ncol(offset) != ncol(Y)) {
        stop("offset must have the same dimensions as Y if it is a matrix.")
      }
    } else if (is.vector(offset) && length(offset) == ncol(Y)) {
      offset <- matrix(rep(offset, n_datasets), nrow = n_datasets, byrow = TRUE)
    } else if (length(offset) == 1) {
      offset <- matrix(offset, nrow = n_datasets, ncol = ncol(Y))
    } else {
      stop("offset must be a matrix, vector, or scalar.")
    }

    # Create processed_data_list
    processed_data_list <- lapply(1:n_datasets, function(i) {
      data.frame(y = Y[i,], x = smooth_var[i,], offset = offset[i,])
    })
  }

  # Validate and process S
  if (!is.null(S)) {
    if (is.list(S)) {
      if (length(S) != n_datasets) {
        stop("S must have the same length as the number of datasets when provided as a list.")
      }
    } else if (is.matrix(S)) {
      S <- split(S, row(S))
    } else if (is.vector(S) & is.numeric(S)) {
      S <- split(matrix(rep(S, n_datasets), nrow = n_datasets, byrow = TRUE), 1:n_datasets)
    } else if(is.character(S)){
      # check if S is a column name in data_list
      if(is.null(data_list)){
        stop("data_list must be provided if S is a character string.")
      }
      if(!all(S %in% names(data_list[[1]]))){
        stop("S must be a valid column name in data_list.")
      }
      S <- lapply(data_list, function(df) df[[S]])
    }
    else {
      stop("S must be a string, list, matrix, or numeric vector.")
    }
  }

  # Validate and process Omega
  if (!is.null(Omega)) {
    if (is.list(Omega)) {
      if (length(Omega) != n_datasets) {
        stop("Omega must have the same length as the number of datasets when provided as a list.")
      }
    } else if (is.matrix(Omega)) {
      Omega <- replicate(n_datasets, Omega, simplify = FALSE)
    } else {
      stop("Omega must be either a list of matrices or a single matrix.")
    }
  }

  # Return structured list
  return(list(data_list = processed_data_list, S = S, Omega = Omega))
}





#' Create a TMB Data Object for a Specific Dataset
#'
#' Processes the `tmbdat` object for a specific dataset, using the provided smoothing variables,
#' standard errors, and precision matrix. This function generates spline basis matrices,
#' penalty matrices, and other quantities required for modeling.
#'
#' @param data_i A single dataset extracted from the `data_list` component of `fash_set_data`.
#'               Must be a list containing `y`, `x`, and `offset`.
#' @param Si A numeric vector representing the standard errors for the dataset. Default is `NULL`.
#' @param Omegai A numeric precision matrix for the dataset. Default is `NULL`.
#' @param num_basis An integer specifying the number of O-Spline basis functions to use for the approximation. Default is 30.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients (`beta`). Default is `1e-6`.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior. Default is 2.
#'
#' @return A list formatted as a `tmbdat` object, containing:
#' \item{y}{A numeric vector of response variables for the dataset.}
#' \item{X}{A sparse matrix representing the design matrix for fixed effects.}
#' \item{P}{A sparse matrix representing the penalty matrix for the spline coefficients.}
#' \item{B}{A sparse matrix representing the design matrix for random effects.}
#' \item{offset}{A numeric vector of offsets for the dataset.}
#' \item{logPdet}{The log determinant of the penalty matrix.}
#' \item{betaprec}{The precision for the fixed effects coefficients.}
#' \item{S}{The standard errors for the dataset, if provided as `Si`.}
#' \item{Omega}{The precision matrix for the dataset as a sparse matrix, if provided as `Omegai`.}
#' \item{log_det_Omega}{The log determinant of the precision matrix, if `Omegai` is provided.}
#'
#' @examples
#' # Example usage
#' Y <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' smooth_var <- matrix(runif(20), nrow = 4, ncol = 5)
#' offset <- 1
#' S <- c(0.5, 0.8, 1.2, 1.0, 0.9)
#' Omega <- diag(5)
#' data <- fash_set_data(Y, smooth_var, offset, S, Omega)
#' tmbdat <- fash_set_tmbdat(data$data_list[[1]], Si = S[[1]], Omegai = Omega)
#'
#' @export
fash_set_tmbdat <- function(data_i, Si = NULL, Omegai = NULL, num_basis = 30, betaprec = 1e-6, order = 2) {
  # Extract smoothing variables and response
  y <- data_i$y
  x <- data_i$x
  offset <- data_i$offset

  # Generate spline knots for the smoothing variable
  knots <- seq(min(x), max(x), length.out = num_basis)

  # Generate the design matrices for the spline basis
  X <- BayesGP:::global_poly_helper(x, p = order)
  P <- BayesGP::compute_weights_precision_helper(knots)
  B <- BayesGP:::local_poly_helper(knots = knots, refined_x = x, p = order)

  # Prepare the tmbdat list
  tmbdat <- list(
    y = y,                                                    # Response variable
    X = as(as(as(X, "dMatrix"), "generalMatrix"), "TsparseMatrix"), # Design matrix for fixed effects
    P = as(as(as(P, "dMatrix"), "generalMatrix"), "TsparseMatrix"), # Penalty matrix for spline coefficients
    B = as(as(as(B, "dMatrix"), "generalMatrix"), "TsparseMatrix"), # Design matrix for random effects
    offset = offset,                                          # Offset term
    logPdet = as.numeric(determinant(P, logarithm = TRUE)$modulus), # Log determinant of penalty matrix
    betaprec = betaprec                                       # Precision for fixed effects
  )

  # Add Omega and compute log determinant if available
  if (!is.null(Si)) {
    # check if Si has the same same length as y
    if (length(Si) != length(y)) {
      # if not, see if Si is a scalar
      if (length(Si) == 1) {
        Si <- rep(Si, length(y))
      } else{
        stop("Standard errors (Si) must have the same length as the response variable (y).")
      }
    }
    tmbdat$S <- Si
  }
  if (!is.null(Omegai)) {
    # check if the dimension of Omegai matches the number of observations
    if (nrow(Omegai) != length(y) || ncol(Omegai) != length(y)) {
      stop("Precision matrix (Omegai) must have the same dimensions as the response variable (y).")
    }
    tmbdat$Omega <- as(as(as(Omegai, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    tmbdat$log_det_Omega <- as.numeric(determinant(Omegai, logarithm = TRUE)$modulus)
  }

  return(tmbdat)
}

