#' Perform Full FASH Analysis
#'
#' This function performs the full FASH pipeline, including data setup, likelihood computation,
#' empirical Bayes estimation, and outputs a structured `fashr` object.
#'
#' @param Y Either a numeric matrix of response variables or a character string specifying the column name in `data_list` for response variables.
#' @param smooth_var A numeric matrix, vector, or a character string specifying the column name in `data_list` for smoothing variables.
#' @param offset A numeric matrix, vector, scalar, or a character string specifying the column name in `data_list` for offset variables. Defaults to 0.
#' @param S A numeric matrix, vector, scalar, or list representing the standard errors of `Y`. Or a character string specifying the column name in `data_list` for SD. Default is `NULL`.
#' @param Omega Either a list of precision matrices (one for each dataset) or a single precision matrix (shared across all datasets). Default is `NULL`.
#' @param data_list A list of data frames, where each data frame corresponds to a single dataset. Default is `NULL`.
#' @param grid A numeric vector representing the grid of PSD (Predictive Standard Deviation) values. Default is `seq(0, 2, length.out = 10)`.
#' @param likelihood A character string specifying the likelihood function to use. Options are `"gaussian"` and `"poisson"`. Default is `"gaussian"`.
#' @param num_basis An integer specifying the number of O-Spline basis functions. Default is 30.
#' @param betaprec A numeric value representing the precision of the fixed effects coefficients (`beta`). Default is `1e-6`.
#' @param order An integer specifying the order of the Integrated Wiener Process (IWP) prior. Default is 2.
#' @param pred_step A numeric value specifying the prediction step size. Default is 1.
#' @param penalty A numeric value representing the lambda value for the Dirichlet prior. Default is 1 (no additional penalty).
#' @param num_cores An integer specifying the number of cores to use for parallel processing. Default is 1.
#' @param verbose A logical value. If `TRUE`, shows progress messages and timing for each step. Default is `FALSE`.
#'
#' @return A `fashr` object containing:
#'   \describe{
#'     \item{\code{prior_weights}}{Estimated prior weights for PSD values.}
#'     \item{\code{posterior_weights}}{Posterior weight matrix of PSD values.}
#'     \item{\code{psd_grid}}{PSD grid values.}
#'     \item{\code{lfdr}}{Local false discovery rate for each dataset.}
#'     \item{\code{settings}}{A list of settings used in the FASH pipeline.}
#'     \item{\code{fash_data}}{A structured list of data components.}
#'   }
#'
#' @examples
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' result <- fashr(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", grid = grid, likelihood = "poisson", verbose = TRUE)
#' print(result)
#'
#' @export
fashr <- function(Y = NULL, smooth_var = NULL, offset = 0, S = NULL, Omega = NULL, data_list = NULL, grid = seq(0, 2, length.out = 10),
                  likelihood = "gaussian", num_basis = 30, betaprec = 1e-6, order = 2, pred_step = 1, penalty = 1,
                  num_cores = 1, verbose = FALSE) {
  # Helper function for timing and verbose output
  timing_message <- function(step_name, code_block) {
    start_time <- Sys.time()
    if (verbose) cat(sprintf("Starting %s...\n", step_name))
    result <- code_block()
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")

    # Format elapsed time dynamically
    if (elapsed_time < 60) {
      time_message <- sprintf("%.2f seconds", elapsed_time)
    } else if (elapsed_time < 3600) {
      time_message <- sprintf("%.2f minutes", elapsed_time / 60)
    } else {
      time_message <- sprintf("%.2f hours", elapsed_time / 3600)
    }

    if (verbose) cat(sprintf("Completed %s in %s.\n", step_name, time_message))
    return(result)
  }

  # Step 1: Data setup
  fash_data <- timing_message("data setup", function() {
    fash_set_data(data_list = data_list, Y = Y, smooth_var = smooth_var, offset = offset, S = S, Omega = Omega)
  })

  # Step 2: Likelihood computation
  L_matrix <- timing_message("likelihood computation", function() {
    fash_L_compute(fash_data, likelihood = likelihood, num_cores = num_cores, grid = grid,
                   num_basis = num_basis, betaprec = betaprec, order = order, pred_step = pred_step,
                   verbose = verbose)
  })

  # Step 3: Empirical Bayes estimation
  eb_result <- timing_message("empirical Bayes estimation", function() {
    fash_eb_est(L_matrix, grid = grid, penalty = penalty)
  })

  # Step 4: Compute additional metrics
  # suppressWarnings is used to avoid warnings from fash_post_ordering
  ordering_lfdr <- suppressWarnings(fash_post_ordering(eb_result, ordering = "lfdr"))

  # Step 5: Create and return the fashr object
  result <- structure(
    list(
      prior_weights = eb_result$prior_weight,
      posterior_weights = eb_result$posterior_weight,
      psd_grid = grid,
      lfdr = (1 - ordering_lfdr$ordered_metrics),
      settings = list(
        num_basis = num_basis,
        betaprec = betaprec,
        order = order,
        pred_step = pred_step,
        likelihood = likelihood,
        penalty = penalty
      ),
      fash_data = fash_data
    ),
    class = "fashr"
  )

  if (verbose) cat("FASHR object created successfully.\n")

  return(result)
}







#' Perform False Discovery Rate (FDR) Control
#'
#' This function performs hypothesis testing by controlling the False Discovery Rate (FDR) based on the
#' local false discovery rate (lfdr) stored in the `fashr` object.
#'
#' @param fashr_obj A `fashr` object containing the results of the FASH pipeline, including `lfdr`.
#' @param alpha A numeric value specifying the significance level for hypothesis testing. Default is `0.05`.
#' @param plot A logical value. If `TRUE`, plots the sorted FDR values with a horizontal line at the `alpha` level. Default is `FALSE`.
#'
#' @return A list containing:
#'
#'       - `fdr_results`: A data frame with columns:
#'
#'       - `index`: The original dataset index.
#'
#'       - `FDR`: The false discovery rate for each dataset.
#'
#' @examples
#' # Example usage
#' data_list <- list(
#'   data.frame(y = rnorm(n = 5, sd = 0.5), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(n = 5, sd = 0.8), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(n = 5, sd = 0.6), x = 1:5, offset = 0),
#'   data.frame(y = rnorm(n = 5, sd = 0.7), x = 1:5, offset = 0)
#' )
#' S <- list(rep(0.5, 5), rep(0.8, 5), rep(0.6, 5), rep(0.7, 5))
#' Omega <- list(diag(1/0.5^2), diag(1/0.8^2), diag(1/0.6^2), diag(1/0.7^2))
#' grid <- seq(0, 2, length.out = 10)
#' fashr_obj <- fashr(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", S = S, Omega = Omega, grid = grid, likelihood = "gaussian", verbose = TRUE)
#' fdr_control(fashr_obj, alpha = 0.05, plot = TRUE)
#'
#' @export
fdr_control <- function(fashr_obj, alpha = 0.05, plot = FALSE) {
  # Extract the local false discovery rates (lfdr)
  lfdr <- fashr_obj$lfdr
  if (is.null(lfdr)) {
    stop("The `fashr` object does not contain local false discovery rates (lfdr).")
  }

  # Sort lfdr and compute the FDR for each dataset
  n <- length(lfdr)
  lfdr_sorted <- sort(lfdr, index.return = TRUE)
  cumulative_lfdr <- cumsum(lfdr_sorted$x) / seq_len(n)

  # Identify significant datasets based on alpha
  significant <- cumulative_lfdr <= alpha
  significant_count <- sum(significant)
  total_count <- n

  # Create results data frame
  fdr_results <- data.frame(
    index = lfdr_sorted$ix,
    FDR = cumulative_lfdr
  )

  # Display message
  message <- sprintf(
    "%d datasets are significant at alpha level %.2f. Total datasets tested: %d.",
    significant_count, alpha, total_count
  )
  cat(message, "\n")

  # Plot the FDR values if plot = TRUE
  if (plot) {
    plot(
      1:n, cumulative_lfdr, type = "b", pch = 19, col = "blue",
      xlab = "Dataset Rank (Sorted by LFDR)", ylab = "Cumulative FDR",
      main = sprintf("FDR Control with Alpha = %.2f", alpha)
    )
    abline(h = alpha, col = "red", lty = 2)  # Horizontal line at alpha
    legend("topright", legend = c("FDR Values", "Alpha Level"), col = c("blue", "red"), lty = c(1, 2), pch = c(19, NA))
  }

  # Return results
  return(list(fdr_results = fdr_results))
}





#' Plot Method for FASHR Objects
#'
#' Generates a structure plot of the posterior weights stored in the `fashr` object,
#' visualizing the distribution of posterior weights across datasets and PSD values.
#'
#' @param x A `fashr` object containing the results of the FASH pipeline.
#' @param ordering A character string specifying the method for reordering datasets (e.g., "mean" or "lfdr").
#'
#'   - `"mean"`: Reorder by the mean of the posterior PSD.
#'
#'   - `"lfdr"`: Reorder by the local false discovery rate.
#'
#'   - `NULL`: No reordering (default).
#'
#' @param discrete A logical value. If `TRUE`, treats PSD values as discrete categories with distinct colors.
#'                 If `FALSE`, treats PSD values as a continuous variable with a gradient. Default is `FALSE`.
#' @param ... Additional arguments passed to `fash_structure_plot`.
#'
#' @examples
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' S <- NULL
#' Omega <- NULL
#' grid <- seq(0, 2, length.out = 10)
#' fashr_obj <- fashr(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", S = S, Omega = Omega, grid = grid, likelihood = "poisson", verbose = TRUE)
#' plot(fashr_obj, ordering = "mean", discrete = TRUE)
#'
#' @export
plot.fashr <- function(x, ordering = NULL, discrete = FALSE, ...) {
  # Validate input
  if (!inherits(x, "fashr")) {
    stop("Input must be a `fashr` object.")
  }

  # Generate the structure plot
  fash_structure_plot(
    eb_output = list(
      posterior_weight = x$posterior_weights,
      prior_weight = x$prior_weights
    ),
    ordering = ordering,
    discrete = discrete,
    ...
  )
}






#' Predict Method for FASHR Objects
#'
#' Generates posterior predictions for a specific dataset from a `fashr` object using Bayesian Model Averaging.
#'
#' @param object A `fashr` object containing the results of the FASH pipeline.
#' @param index An integer specifying the dataset index to predict. Default is `1`.
#' @param smooth_var A numeric vector specifying refined x values for prediction. If `NULL`, uses the x values from the model fit. Default is `NULL`.
#' @param only.samples A logical value. If `TRUE`, returns posterior samples. If `FALSE`, summarizes the samples into mean and 95 percent confidence intervals. Default is `FALSE`.
#' @param M An integer specifying the number of posterior samples to generate. Default is `3000`.
#' @param ... Additional arguments (not used).
#'
#' @return If `only.samples = TRUE`, a matrix of posterior samples where rows correspond to `smooth_var` and columns correspond to posterior draws.
#' If `only.samples = FALSE`, a data frame summarizing posterior predictions with columns:
#'
#'   - `x`: The refined x values.
#'
#'   - `mean`: The posterior mean.
#'
#'   - `lower`: The lower bound of the 95 percent interval.
#'
#'   - `upper`: The upper bound of the 95 percent interval.
#'
#'   - `median`: The posterior median.
#'
#' @examples
#' # Example 1: Predict for a specific dataset with summarized results
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' S <- list(rep(0.5, 5), rep(0.8, 5))
#' Omega <- list(diag(5), diag(5))
#' grid <- seq(0, 2, length.out = 10)
#' fashr_obj <- fashr(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", S = S, Omega = Omega, grid = grid, likelihood = "poisson", verbose = TRUE)
#' predict(fashr_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = FALSE)
#'
#' # Example 2: Generate posterior samples
#' samples <- predict(fashr_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = TRUE)
#' dim(samples)  # Rows: refined_x, Columns: posterior samples
#'
#' # Example 3: Use original x values for prediction
#' summary <- predict(fashr_obj, index = 1, smooth_var = NULL, only.samples = FALSE)
#' head(summary)
#'
#' # Example 4: Increase number of posterior samples
#' samples <- predict(fashr_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = TRUE, M = 5000)
#' summary <- predict(fashr_obj, index = 1, smooth_var = seq(1, 5, length.out = 50), only.samples = FALSE, M = 5000)
#'
#' @export
predict.fashr <- function(object, index = 1, smooth_var = NULL, only.samples = FALSE, M = 3000, ...) {
  # Validate input
  if (!inherits(object, "fashr")) {
    stop("Input must be a `fashr` object.")
  }
  if (index < 1 || index > length(object$posterior_weights)) {
    stop("Index is out of range for the datasets in the `fashr` object.")
  }

  # Extract dataset-specific components
  data_i <- object$fash_data$data_list[[index]]
  Si <- object$fash_data$S[[index]]
  Omegai <- object$fash_data$Omega[[index]]
  psd_values <- object$prior_weights$psd
  posterior_weights <- object$posterior_weights[index, ]

  # Use smooth_var if provided; otherwise, default to the dataset's x values
  refined_x <- if (!is.null(smooth_var)) smooth_var else data_i$x

  # Retrieve settings from fashr object
  settings <- object$settings

  # Generate posterior samples using BMA
  posterior_samples <- fash_bma_sampling(
    data_i = data_i,
    posterior_weights = posterior_weights,
    psd_values = psd_values,
    refined_x = refined_x,
    M = M,
    Si = Si,
    Omegai = Omegai,
    num_basis = settings$num_basis,
    betaprec = settings$betaprec,
    order = settings$order,
    pred_step = settings$pred_step,
    likelihood = settings$likelihood
  )$posterior_samples

  # Return samples if only.samples = TRUE
  if (only.samples) {
    return(posterior_samples)
  }

  # Summarize samples: mean and 95% intervals
  posterior_summary <- data.frame(
    x = refined_x,
    mean = rowMeans(posterior_samples),
    median = apply(posterior_samples, 1, median),
    lower = apply(posterior_samples, 1, function(x) quantile(x, 0.025)),
    upper = apply(posterior_samples, 1, function(x) quantile(x, 0.975))
  )

  return(posterior_summary)
}





#' Print Method for FASHR Objects
#'
#' Displays a summary of the fitted `fashr` object, including the number of datasets,
#' type of likelihood used, the number of PSD grid values, and the order of the Integrated Wiener Process (IWP).
#'
#' @param x A `fashr` object.
#' @param ... Additional arguments (not used).
#'
#' @examples
#' data_list <- list(
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0),
#'   data.frame(y = rpois(5, lambda = 5), x = 1:5, offset = 0)
#' )
#' grid <- seq(0, 2, length.out = 10)
#' fashr_obj <- fashr(data_list = data_list, Y = "y", smooth_var = "x", offset = "offset", grid = grid, likelihood = "poisson", verbose = TRUE)
#' print(fashr_obj)
#'
#' @export
print.fashr <- function(x, ...) {
  # Validate input
  if (!inherits(x, "fashr")) {
    stop("Input must be a `fashr` object.")
  }

  # Extract relevant information
  n_datasets <- length(x$fash_data$data_list)
  likelihood <- x$settings$likelihood
  n_grid_initial <- length(x$psd_grid)
  nontrivial_grid_values <- sum(x$prior_weights$prior_weight > 0)
  iwp_order <- x$settings$order

  # Display summary
  cat("Fitted FASHR Object\n")
  cat("-------------------\n")
  cat(sprintf("Number of datasets: %d\n", n_datasets))
  cat(sprintf("Likelihood: %s\n", likelihood))
  cat(sprintf("Number of PSD grid values: %d (initial), %d (non-trivial)\n", n_grid_initial, nontrivial_grid_values))
  cat(sprintf("Order of Integrated Wiener Process (IWP): %d\n", iwp_order))
}






