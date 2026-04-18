# Quick health check on the results table
raw <- read.csv("benchmark_results/raw_results.csv", stringsAsFactors = FALSE)

cat(sprintf("Total rows:      %d\n", nrow(raw)))
cat(sprintf("Expected:        %d  (8 cells x 6 models x 10 folds)\n",
            8 * 6 * 10))

cat("\n--- Rows per (dataset, regime, model) ---\n")
tab <- table(raw$dataset, raw$regime, raw$model)
# Flatten to a readable format
for (ds in dimnames(tab)[[1]]) {
  for (reg in dimnames(tab)[[2]]) {
    for (m in dimnames(tab)[[3]]) {
      n <- tab[ds, reg, m]
      if (n > 0 && n != 10) {
        cat(sprintf("  !! %s @ %s :: %s  n=%d (expected 10)\n",
                    ds, reg, m, n))
      }
    }
  }
}
cat("(any lines printed above indicate folds that did not complete)\n")

cat("\n--- Fit status breakdown ---\n")
print(table(raw$dataset, raw$fit_status))

cat("\n--- Non-convergence flags ---\n")
print(table(raw$dataset, raw$model, raw$converged, useNA = "ifany"))

cat("\n--- NA rates per metric, by dataset ---\n")
for (metric in c("harrell_c", "uno_c", "ibs", "calib_slope")) {
  na_by_ds <- tapply(is.na(raw[[metric]]), raw$dataset, sum)
  cat(sprintf("  %s: %s\n", metric,
              paste(names(na_by_ds), na_by_ds, sep = "=", collapse = ", ")))
}



#----------------------------------------------------------------------

agg <- read.csv("benchmark_results/aggregated.csv", stringsAsFactors = FALSE)

# Harrell C by model, averaged across cells
cat("=== Median Harrell C per (dataset, regime, model) ===\n")
print(
  agg[order(agg$dataset, agg$regime, -agg$harrell_c_median),
      c("dataset", "regime", "model",
        "harrell_c_median", "harrell_c_q25", "harrell_c_q75",
        "harrell_c_n")],
  row.names = FALSE
)

cat("\n=== Median IBS per (dataset, regime, model) ===\n")
print(
  agg[order(agg$dataset, agg$regime, agg$ibs_median),
      c("dataset", "regime", "model",
        "ibs_median", "ibs_q25", "ibs_q75")],
  row.names = FALSE
)

cat("\n=== Partition (best model per cell by Harrell C) ===\n")
part <- read.csv("benchmark_results/partition.csv", stringsAsFactors = FALSE)
print(part[, c("dataset", "regime", "c_observed", "epv_observed",
               "best_model", "best_metric_median", "gap_to_second")],
      row.names = FALSE)




#-----------------------------------------------------------

res <- read.csv("benchmark_results/raw_results.csv", stringsAsFactors = FALSE)

# Per-cell medians for metrics not shown yet
library(dplyr)
cell_summary <- res %>%
  group_by(dataset, regime, model) %>%
  summarise(
    calib_slope_med = median(calib_slope, na.rm = TRUE),
    uno_c_med       = median(uno_c, na.rm = TRUE),
    fit_time_med    = median(fit_time_sec, na.rm = TRUE),
    n_folds         = sum(is.finite(harrell_c)),
    .groups = "drop"
  )
print(cell_summary, n = Inf)

# Pooled medians across ALL cells, for Table I
overall <- res %>%
  group_by(model) %>%
  summarise(
    harrell_c   = median(harrell_c,   na.rm = TRUE),
    uno_c       = median(uno_c,       na.rm = TRUE),
    ibs         = median(ibs,         na.rm = TRUE),
    calib_slope = median(calib_slope, na.rm = TRUE),
    fit_time    = median(fit_time_sec, na.rm = TRUE)
  )
print(overall)

# Pairwise significance (if file exists)
if (file.exists("benchmark_results/pairwise_wilcoxon.csv")) {
  pw <- read.csv("benchmark_results/pairwise_wilcoxon.csv")
  print(table(pw$metric, pw$conclusion))
}

# PH test results
ph <- res %>% group_by(dataset, regime) %>%
  summarise(ph_p_med = median(ph_global_p, na.rm = TRUE), .groups = "drop")
print(ph)
