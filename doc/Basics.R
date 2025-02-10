## ----setup--------------------------------------------------------------------
knitr::opts_chunk$set(fig.width = 8, fig.height = 6)
library(fashr)

## ----echo=FALSE---------------------------------------------------------------
### Write a function, to simulate a dataset with size n, given a function g
simulate_data <- function(g, sd=0.1){
  x <- 1:16
  # simulate sd from sampling from sd with replacement
  sd <- sample(x = sd, size = length(x), replace = TRUE)
  y <- g(x) + rnorm(n = length(x), sd = sd)
  return(data.frame(x=x, y=y,truef = g(x), sd = sd))
}

### Write a function, to simulate a random function as g using cubic B-spline basis representation with random basis weights:
simulate_nonlinear_function <- function(n_basis = 20, sd_function = 1, sd_poly = 0.1) {
  if(n_basis < 3) stop("n_basis must be greater than 3")
  # Define the range and knots for the B-spline basis
  x_min <- 0
  x_max <- 16

  # Generate equally spaced knots within the range [x_min, x_max]
  knots <- seq(x_min, x_max, length.out = n_basis - 3)

  # Generate random weights for the basis functions
  pred_step = 1
  p = 1
  sd_function <- sd_function/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))
  prec_mat <- (1/sd_function^2) * BayesGP:::compute_weights_precision_helper(knots)
  weights <- as.vector(LaplacesDemon::rmvnp(n = 1, mu = rep(0, ncol(prec_mat)), Omega = prec_mat))
  # Generate random weights for the linear functions
  beta_vec <- rnorm(n = p, mean = 0, sd = sd_poly)

  # Return a function that evaluates the spline at new x values
  function(x_new) {
    # Create the B-spline basis for the new x values using the predefined knots
    spline_new <- BayesGP:::local_poly_helper(knots = knots, refined_x = x_new, p = p)
    x_new_design <- BayesGP:::global_poly_helper(x = x_new, p = p)
    # Return the function
    return(x_new_design %*% beta_vec + as.vector(spline_new %*% weights))

  }
}

### Write a function, to simulate a random function as g using constant and linear function
simulate_linear_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  # beta1 <- rnorm(1, mean = 0, sd = sd_poly)
  beta1 <- sample(c(-0.5,0.5),1)
  function(x_new) {
    return(beta0 + beta1 * x_new)
  }
}

### Write a function, to simulate a random function as g using constant, linear and quadratic function
simulate_quadratic_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  beta1 <- rnorm(1, mean = 0, sd = sd_poly)
  beta2 <- rnorm(1, mean = 0, sd = sd_poly)
  function(x_new) {
    return(beta0 + beta1 * x_new + beta2 * x_new^2)
  }
}

### Write a function, to simulate a nondynamic function as g using constant function
simulate_nondynamic_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  function(x_new) {
    return(beta0)
  }
}

### simulate process: first draw a random function g, then draw a dataset from g
simulate_process <- function(n_basis = 50, sd_fun = 1, sd = 0.1, sd_poly = 0.1, type = "nonlinear"){
  if(type == "linear"){
    g <- simulate_linear_function(sd_poly = sd_poly)

  }
  else if(type == "quadratic"){
    g <- simulate_quadratic_function(sd_poly = sd_poly)
  }
  else if(type == "nonlinear") {
    g <- simulate_nonlinear_function(n_basis = n_basis, sd_function = sd_fun, sd_poly = sd_poly)
  }
  else if(type == "nondynamic") {
    g <- simulate_nondynamic_function(sd_poly = sd_poly)
  }
  else {
    stop("type must be one of 'linear', 'nonlinear', 'nondynamic'")
  }
  
  return(simulate_data(g = g, sd = sd))
}


## ----echo=FALSE---------------------------------------------------------------
num_cores <- 4
set.seed(1234)
N <- 100
propA <- 0.5; propB <- 0.3; propC <- 0.2
sigma_vec <- c(0.1, 0.3, 0.5, 1)

sizeA <- N * propA
data_sim_list_A <- lapply(1:sizeA, function(i) simulate_process(sd_poly = 1, type = "nondynamic", sd = sigma_vec))

sizeB <- N * propB
if(sizeB > 0){
data_sim_list_B <- lapply(1:sizeB, function(i) simulate_process(sd_poly = 0.5, type = "linear", sd = sigma_vec))

}else{
  data_sim_list_B <- list()
}

sizeC <- N * propC
data_sim_list_C <- lapply(1:sizeC, function(i) simulate_process(sd_poly = 0.1, type = "nonlinear", sd = sigma_vec, sd_fun = 1))

datasets <- c(data_sim_list_A, data_sim_list_B, data_sim_list_C)
sigma <- unlist(lapply(datasets, function(x) unique(x$sd)))

labels <- c(rep("A", sizeA), rep("B", sizeB), rep("C", sizeC))

## -----------------------------------------------------------------------------
length(datasets)
str(datasets[[1]])

## -----------------------------------------------------------------------------
table(labels)

## -----------------------------------------------------------------------------
fash_fit <- fashr(Y = "y", smooth_var = "x", S = "sd", data_list = datasets, 
                  likelihood = "gaussian", order = 2,
                  num_cores = 2, num_basis = 20, grid = seq(0, 2, length.out = 20),
                  verbose = TRUE)

## -----------------------------------------------------------------------------
fash_fit

## -----------------------------------------------------------------------------
plot(fash_fit)

## -----------------------------------------------------------------------------
plot(fash_fit, discrete = TRUE, ordering = "lfdr")

## -----------------------------------------------------------------------------
fdr_result <- fdr_control(fash_fit, alpha = 0.1, plot = TRUE)

## -----------------------------------------------------------------------------
detected_indices <- fdr_result$fdr_results$index[1:22]

## -----------------------------------------------------------------------------
sum(labels[detected_indices] == "C")/20

## -----------------------------------------------------------------------------
sum(labels[detected_indices] != "C")/22

## -----------------------------------------------------------------------------
fitted_beta <- predict(fash_fit, index = detected_indices[1])
str(fitted_beta)

## -----------------------------------------------------------------------------
fitted_beta_new <- predict(fash_fit, index = detected_indices[1], smooth_var = seq(0, 16, length.out = 100))
str(fitted_beta_new)

## -----------------------------------------------------------------------------
fitted_beta_samples <- predict(fash_fit, index = detected_indices[1], 
                               smooth_var = seq(0, 16, length.out = 100), 
                               only.samples = TRUE, M = 50)
str(fitted_beta_samples)

## -----------------------------------------------------------------------------
plot(datasets[[detected_indices[1]]]$x, datasets[[detected_indices[1]]]$y, type = "p", col = "black", lwd = 2, xlab = "Time", ylab = "Effect Size")
lines(fitted_beta_new$x, fitted_beta_new$mean, col = "red", lwd = 2)
polygon(c(fitted_beta_new$x, rev(fitted_beta_new$x)), c(fitted_beta_new$lower, rev(fitted_beta_new$upper)), col = rgb(1, 0, 0, 0.2), border = NA)

## -----------------------------------------------------------------------------
fash_fit_2 <- fashr(Y = "y", smooth_var = "x", S = "sd", data_list = datasets, 
                  likelihood = "gaussian", order = 1,
                  num_cores = 2, num_basis = 20, grid = seq(0, 2, length.out = 20),
                  verbose = TRUE)

## -----------------------------------------------------------------------------
fash_structure_plot(fash_fit_2, discrete = TRUE, ordering = "lfdr")

## -----------------------------------------------------------------------------
fdr_result_2 <- fdr_control(fash_fit_2, alpha = 0.1, plot = TRUE)
detected_indices_2 <- fdr_result_2$fdr_results$index[1:55]

## -----------------------------------------------------------------------------
sum(labels[detected_indices_2] != "A")/50

## -----------------------------------------------------------------------------
sum(labels[detected_indices_2] == "A")/55

