#' Perform Empirical Bayes (EB) on Inferred Likelihood Matrix
#'
#' This function performs empirical Bayes estimation on the provided (log) likelihood matrix, incorporating an optional Dirichlet prior penalty.
#'
#' @param L_matrix A numeric matrix representing the log-likelihood values. Rows correspond to datasets, and columns correspond to grid points.
#' @param penalty A numeric value representing the lambda value for the Dirichlet prior.
#' @param grid A numeric vector representing the grid of PSD (Predictive Standard Deviation) values.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{prior_weight}}{A data frame with two columns: \code{psd} (non-trivial PSD values) and \code{prior_weight} (learned prior weights).}
#'     \item{\code{posterior_weight}}{A numeric matrix where rows correspond to datasets, and columns correspond to non-trivial PSD values, representing the posterior weights.}
#'   }
#'
#' @examples
#' # Example usage
#' L_matrix <- matrix(rnorm(500), nrow = 100, ncol = 5)
#' grid <- seq(0, 2, length.out = 5)
#' result <- fash_eb_est(L_matrix, penalty = 2, grid = grid)
#' print(result$prior_weight)
#' print(result$posterior_weight)
#'
#' @importFrom mixsqp mixsqp
#'
#' @export
#'
#'
fash_eb_est <- function(L_matrix, penalty = 1, grid) {
  num_datasets <- nrow(L_matrix)
  num_components <- ncol(L_matrix)

  # Add a Dirichlet prior penalty if specified
  if (penalty > 1) {
    prior_null <- matrix(0, nrow = penalty - 1, ncol = num_components)
    prior_null[, 1] <- 1  # Prior mass on the first grid point
    L_matrix_original <- rbind(exp(L_matrix), prior_null)
    fit.sqp <- mixsqp::mixsqp(
      L = L_matrix_original,
      log = FALSE,
      control = list(tol.svd = 0, verbose = FALSE)
    )
  } else {
    fit.sqp <- mixsqp::mixsqp(
      L = L_matrix,
      log = TRUE,
      control = list(tol.svd = 0, verbose = FALSE)
    )
  }

  # Extract non-trivial PSD values and corresponding prior weights
  non_trivial <- which(fit.sqp$x > 0)
  prior_weight <- data.frame(
    psd = grid[non_trivial],
    prior_weight = fit.sqp$x[non_trivial]
  )

  # Compute posterior weights for each dataset
  posterior_weight <- matrix(0, nrow = num_datasets, ncol = length(non_trivial))
  for (i in 1:num_datasets) {
    exp_values <- exp(L_matrix[i, ] - max(L_matrix[i, ]) + log(fit.sqp$x))
    normalized_values <- exp_values[non_trivial] / sum(exp_values[non_trivial])
    posterior_weight[i, ] <- normalized_values
  }
  colnames(posterior_weight) <- as.character(grid[non_trivial])

  # Return results
  return(list(
    prior_weight = prior_weight,
    posterior_weight = posterior_weight
  ))
}




#' Structure Plot for Posterior Weights
#'
#' This function takes the output of \code{fash_eb_est} and generates a structure plot
#' visualizing the posterior weights for all datasets. It can display PSD values
#' as either continuous or discrete variables and optionally reorder datasets.
#'
#' @param eb_output A list output from \code{fash_eb_est}, containing:
#'   \describe{
#'     \item{posterior_weight}{A numeric matrix of posterior weights (datasets as rows, PSD as columns).}
#'     \item{prior_weight}{A data frame of prior weights (not used in this plot).}
#'   }
#' @param discrete A logical value. If \code{TRUE}, treats PSD values as discrete categories with distinct colors.
#'                 If \code{FALSE}, treats PSD values as a continuous variable with a gradient.
#' @param ordering A character string specifying the method for reordering datasets. Options are:
#'   \describe{
#'     \item{NULL}{No reordering (default).}
#'     \item{`mean`}{Reorder by the mean of the posterior PSD.}
#'     \item{`median`}{Reorder by the median of the posterior PSD.}
#'     \item{`lfdr`}{Reorder by the local false discovery rate (posterior probability of PSD = 0).}
#'   }
#' @param selected_indices A numeric vector specifying the indices of datasets to display. If \code{NULL}, displays all datasets.
#' @return A ggplot object representing the structure plot.
#'
#' @examples
#' # Example usage
#' grid <- seq(0.1, 2, length.out = 5)
#' L_matrix <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' eb_output <- fash_eb_est(L_matrix, penalty = 2, grid = grid)
#' plot_cont <- fash_structure_plot(eb_output, discrete = FALSE, ordering = "mean")
#' plot_disc <- fash_structure_plot(eb_output, discrete = TRUE, ordering = "median")
#' print(plot_cont)
#' print(plot_disc)
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs scale_fill_brewer scale_fill_gradient coord_flip theme_minimal theme element_blank element_rect
#' @importFrom reshape2 melt
#'
#' @export
fash_structure_plot <- function(eb_output, discrete = FALSE, ordering = NULL, selected_indices = NULL) {

  # Select indices if specified
  if (!is.null(selected_indices)) {
    eb_output_selected <- eb_output
    eb_output_selected$posterior_weight <- eb_output$posterior_weight[selected_indices, , drop = FALSE]
  } else {
    eb_output_selected <- eb_output
  }

  # Extract posterior weights matrix
  posterior_weights_matrix <- eb_output$posterior_weight

  # Reorder datasets if ordering is specified
  if (!is.null(ordering)) {
    order_result <- fash_post_ordering(eb_output_selected, ordering = ordering)
    posterior_weights_matrix <- order_result$ordered_matrix
    ordered_indices <- order_result$ordered_indices
  } else {
    ordered_indices <- seq_len(nrow(posterior_weights_matrix))
  }

  # Extract PSD values
  psd_values <- as.numeric(colnames(posterior_weights_matrix))

  # Convert the posterior matrix to a data frame for ggplot
  posterior_weights_df <- as.data.frame(posterior_weights_matrix)
  posterior_weights_df$id <- ordered_indices

  # Melt the data frame for ggplot
  melted_data <- reshape2::melt(posterior_weights_df, id.vars = "id")
  melted_data$variable <- as.numeric(as.character(melted_data$variable))

  # Adjust the PSD variable for discrete or continuous plotting
  if (discrete) {
    # Round PSD values and convert to factor
    melted_data$variable <- factor(round(melted_data$variable, 3), levels = round(psd_values, 3))
    fill_scale <- ggplot2::scale_fill_brewer(palette = "Set3", name = "PSD (Rounded)")
  } else {
    fill_scale <- ggplot2::scale_fill_gradient(low = "white", high = "blue", name = "PSD")
  }

  # Create the structure plot
  ggplot2::ggplot(melted_data, ggplot2::aes(x = factor(id, levels = posterior_weights_df$id), y = value, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(
      x = "Datasets",
      y = "Posterior Weight",
      title = "Structure Plot of Posterior Weights"
    ) +
    fill_scale +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background = ggplot2::element_rect(fill = "white")
    )
}







#' Order Posterior Weight Matrix
#'
#' Orders the posterior weight matrix from \code{fash_eb_est} based on a specified method of the posterior PSD.
#'
#' @param eb_output A list output from \code{fash_eb_est}, containing:
#'   \describe{
#'     \item{\code{posterior_weight}}{A numeric matrix of posterior weights (datasets as rows, PSD as columns).}
#'     \item{\code{prior_weight}}{A data frame of prior weights (not used in this function).}
#'   }
#' @param ordering A character string specifying the method for reordering datasets. Options are:
#'   \describe{
#'     \item{\code{"mean"}}{Reorder by the mean of the posterior PSD (default).}
#'     \item{\code{"median"}}{Reorder by the median of the posterior PSD.}
#'     \item{\code{"lfdr"}}{Reorder by the local false discovery rate (posterior probability of PSD = 0).}
#'   }
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{ordered_matrix}}{The reordered posterior weight matrix.}
#'     \item{\code{ordered_indices}}{The indices used to reorder the matrix.}
#'     \item{\code{ordered_metrics}}{The metrics used for ordering, aligned with the reordered matrix.}
#'   }
#'
#' @examples
#' # Example usage
#' grid <- seq(0.1, 2, length.out = 5)
#' L_matrix <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' eb_output <- fash_eb_est(L_matrix, penalty = 2, grid = grid)
#' result <- fash_post_ordering(eb_output, ordering = "mean")
#' print(result$ordered_matrix)
#' print(result$ordered_indices)
#' print(result$ordered_metrics)
#'
#' @export
fash_post_ordering <- function(eb_output, ordering = "mean") {
  # Extract posterior weight matrix and PSD values
  posterior_weight <- eb_output$posterior_weight
  psd_values <- as.numeric(colnames(posterior_weight))

  # Check if PSD = 0 column exists
  if (ordering == "lfdr" && !0 %in% psd_values) {
    warning("PSD = 0 is not present in the posterior matrix. Returning original ordering.")
    return(list(
      ordered_matrix = posterior_weight,
      ordered_indices = seq_len(nrow(posterior_weight)),
      ordered_metrics = rep(0, nrow(posterior_weight))
    ))
  }

  # Compute ordering metric
  ordering_metric <- switch(ordering,
                            "mean" = apply(posterior_weight, 1, function(row) sum(row * psd_values)),
                            "median" = apply(posterior_weight, 1, function(row) {
                              cum_weights <- cumsum(row)
                              psd_values[which(cum_weights >= 0.5)[1]]
                            }),
                            "lfdr" = (1 - posterior_weight[, which(psd_values == 0)]),
                            stop("Invalid ordering. Choose from 'mean', 'median', or 'lfdr'.")
  )

  # Order rows based on the computed metric
  ordered_indices <- order(ordering_metric, decreasing = FALSE)
  ordered_matrix <- posterior_weight[ordered_indices, , drop = FALSE]
  ordered_metrics <- ordering_metric[ordered_indices]

  # Return the ordered matrix, indices, and metrics
  return(list(
    ordered_matrix = ordered_matrix,
    ordered_indices = ordered_indices,
    ordered_metrics = ordered_metrics
  ))
}



