# Load required libraries
library(microbenchmark)

# Generate synthetic data for testing
set.seed(123)
n_datasets <- 100  # Adjust for larger testing
n_points <- 30

# Create example datasets
data_list <- lapply(1:n_datasets, function(i) {
  df <- data.frame(x = seq(from = 0, to = 3, length.out = n_points), offset = 0)
  df$y <- rpois(n_points, lambda = exp(3 + rnorm(1, sd = 0.1) * df$x))
  df
})

# Define PSD grid values
grid <- seq(0, 0.1, length.out = 10)

# Fit the fash object
fash_obj <- fash(data_list = data_list, Y = "y", smooth_var = "x", grid = grid, likelihood = "poisson", verbose = TRUE, order = 1, penalty = 2)

# Set the number of cores for parallel execution
num_cores <- 1  # Adjust based on your machine's capabilities

# Benchmarking min_lfsr_sampling (sampling-based method)
runtime_sampling <- microbenchmark(
  min_lfsr_sampling(fash_obj, num_cores = num_cores, M = 3000),
  times = 5  # Number of repetitions for averaging runtime
)

# Benchmarking min_lfsr_summary (analytical method)
runtime_summary <- microbenchmark(
  min_lfsr_summary(fash_obj, num_cores = num_cores),
  times = 5  # Number of repetitions for averaging runtime
)

# Print runtime comparisons
print(runtime_sampling)
print(runtime_summary)

# Compare execution times visually
boxplot(runtime_sampling$time / 1e9, runtime_summary$time / 1e9, 
        names = c("Sampling", "Summary"), 
        ylab = "Time (seconds)", main = "Runtime Comparison: Sampling vs Summary")
