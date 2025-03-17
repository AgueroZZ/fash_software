library(fashr)

load("~/Desktop/FASH/FASHresultsummary/output/dynamic_eQTL_real/fash_fit1_all.RData")
load("~/Desktop/FASH/FASHresultsummary/output/dynamic_eQTL_real/fash_fit2_all.RData")

datasets <- fash_fit1$fash_data$data_list
for (i in 1:length(datasets)) {
  datasets[[i]]$SE <- fash_fit1$fash_data$S[[i]]
}
all_genes <- unique(sapply(strsplit(names(datasets), "_"), "[[", 1))



load("~/Desktop/FASH/FASHresultsummary/output/dynamic_eQTL_real/min_lfsr_summary1.RData")
min_lfsr_summary1$index[1:5]



selected_index <- 250533
pred_result <- predict(fash_fit2, index = selected_index, deriv = 0)
plot(pred_result$median ~ pred_result$x, type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "y",
     ylim = range(c(pred_result$lower-0.5, pred_result$upper+0.5)))
polygon(c(pred_result$x, rev(pred_result$x)), c(pred_result$lower, rev(pred_result$upper)), col = rgb(0, 0, 1, 0.3), border = NA)
points(datasets[[selected_index]]$x, datasets[[selected_index]]$y, col = "red", pch = 19)
arrows(
  datasets[[selected_index]]$x,
  datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE,
  datasets[[selected_index]]$x,
  datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE,
  length = 0.05,
  angle = 90,
  code = 3,
  col = "red"
)


lfsr_summary <- compute_lfsr_summary(fash_fit1, index = selected_index, deriv = 0)
min_lfsr_summary <- min(lfsr_summary$lfsr)
min_lfsr_summary
lfsr_sampling <- compute_lfsr_sampling(fash_fit1, index = selected_index, deriv = 0)
min_lfsr_sampling <- min(lfsr_sampling$lfsr)
min_lfsr_sampling





