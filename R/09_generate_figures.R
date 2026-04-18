## ============================================================
## generate_figures.R
##
## Generates all figures for the benchmark paper, reading from
## benchmark_results/raw_results.csv. Each figure is:
##   - displayed in the R graphics device
##   - saved as a high-resolution PNG in benchmark_results/figures/
##
## Figures produced:
##   fig1_c_index.png        C-index vs censoring rate (3 panels)
##   fig2_ibs.png            IBS vs censoring rate (3 panels)
##   fig3_calib_slope.png    Calibration slope grouped bars (3 panels)
##   fig4_partition.png      Empirical partition of (c, EPV) plane
##   fig5_fit_time.png       Wall-clock training time comparison
##
## Requirements: base R only (no ggplot2).
## ============================================================

## ------------------------------------------------------------
## 0. Setup
## ------------------------------------------------------------
res <- read.csv("benchmark_results/raw_results.csv",
                stringsAsFactors = FALSE)

fig_dir <- "benchmark_results/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## --- Model display order, labels, and visual encoding ---
model_order  <- c("cox", "lasso_cox", "ben_cox", "rsf", "deepsurv", "cox_time")
model_labels <- c("Cox", "Lasso-Cox", "BEN-Cox", "RSF", "DeepSurv", "Cox-Time")

model_colors <- c(
  cox       = "#888888",
  lasso_cox = "#2166ac",
  ben_cox   = "#b2182b",
  rsf       = "#1b7837",
  deepsurv  = "#e08214",
  cox_time  = "#c51b7d"
)
model_pch <- c(
  cox       = 15,   # filled square
  lasso_cox = 18,   # filled diamond
  ben_cox   = 16,   # filled circle
  rsf       = 17,   # filled triangle up
  deepsurv  = 25,   # inverted triangle
  cox_time  = 8     # star/asterisk
)
model_lty <- c(
  cox       = 2,    # dashed
  lasso_cox = 1,    # solid
  ben_cox   = 1,
  rsf       = 1,
  deepsurv  = 2,
  cox_time  = 2
)

## --- Dataset display names and censoring levels ---
## Build the mapping from (dataset, regime) to observed c
cell_info <- aggregate(
  cbind(c_achieved, epv, epv_raw) ~ dataset + regime,
  data = res,
  FUN  = function(x) median(x, na.rm = TRUE)
)

## Canonical ordering of datasets for panels
ds_order <- c("SUPPORT", "METABRIC", "TCGA-BRCA")

## --- Aggregate medians and IQR per (dataset, regime, model) ---
agg <- do.call(rbind, lapply(
  split(res, interaction(res$dataset, res$regime, res$model, drop = TRUE)),
  function(g) {
    data.frame(
      dataset = g$dataset[1],
      regime  = g$regime[1],
      model   = g$model[1],
      c_obs   = median(g$c_achieved, na.rm = TRUE),
      epv_obs = median(g$epv, na.rm = TRUE),
      harrell_c_med = median(g$harrell_c, na.rm = TRUE),
      harrell_c_q25 = quantile(g$harrell_c, 0.25, na.rm = TRUE),
      harrell_c_q75 = quantile(g$harrell_c, 0.75, na.rm = TRUE),
      ibs_med       = median(g$ibs, na.rm = TRUE),
      ibs_q25       = quantile(g$ibs, 0.25, na.rm = TRUE),
      ibs_q75       = quantile(g$ibs, 0.75, na.rm = TRUE),
      calib_med     = median(g$calib_slope, na.rm = TRUE),
      uno_c_med     = median(g$uno_c, na.rm = TRUE),
      fit_time_med  = median(g$fit_time_sec, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
rownames(agg) <- NULL


## ------------------------------------------------------------
## Helper: save current plot as PNG and display it
## ------------------------------------------------------------
save_and_show <- function(filename, width = 10, height = 4, res = 300) {
  ## Record the current plot, save to PNG, then replay on screen
  p <- recordPlot()
  png(file.path(fig_dir, filename), width = width, height = height,
      units = "in", res = res)
  replayPlot(p)
  dev.off()
  cat(sprintf("Saved: %s/%s\n", fig_dir, filename))
}


## ------------------------------------------------------------
## Helper: three-panel line plot with error bars
## ------------------------------------------------------------
plot_metric_panels <- function(metric_med, metric_q25, metric_q75,
                               ylab, main_prefix,
                               higher_better = TRUE) {
  ## Layout: 3 panels side by side + shared legend at bottom
  layout_mat <- matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE)
  layout(layout_mat, heights = c(5, 1))
  par(mar = c(4, 4.2, 2.5, 0.8), family = "serif", cex = 0.85)
  
  for (ds in ds_order) {
    sub <- agg[agg$dataset == ds, ]
    ## Sort by observed censoring rate
    c_vals <- sort(unique(sub$c_obs))
    n_c <- length(c_vals)
    
    ## Y limits across all models in this panel
    ymin <- min(sub[[metric_q25]], na.rm = TRUE)
    ymax <- max(sub[[metric_q75]], na.rm = TRUE)
    ypad <- (ymax - ymin) * 0.15
    ylim <- c(ymin - ypad, ymax + ypad)
    
    ## X limits with padding
    if (n_c > 1) {
      xpad <- (max(c_vals) - min(c_vals)) * 0.18
      xlim <- c(min(c_vals) - xpad, max(c_vals) + xpad)
    } else {
      xlim <- c(c_vals[1] - 0.06, c_vals[1] + 0.06)
    }
    
    plot(NA, xlim = xlim, ylim = ylim,
         xlab = expression("Censoring rate " * italic(c)),
         ylab = if (ds == ds_order[1]) ylab else "",
         main = paste0("(", letters[which(ds_order == ds)], ") ", ds),
         xaxt = "n", las = 1)
    axis(1, at = c_vals, labels = sprintf("%.2f", c_vals))
    grid(col = "#e0e0e0", lty = 1, lwd = 0.5)
    
    ## Slight horizontal jitter so error bars don't overlap
    jitter_offsets <- seq(-0.008, 0.008, length.out = length(model_order))
    
    for (i in seq_along(model_order)) {
      m <- model_order[i]
      ms <- sub[sub$model == m, ]
      ms <- ms[order(ms$c_obs), ]
      x <- ms$c_obs + jitter_offsets[i]
      y <- ms[[metric_med]]
      lo <- ms[[metric_q25]]
      hi <- ms[[metric_q75]]
      
      if (length(x) > 1) {
        lines(x, y, col = model_colors[m], lty = model_lty[m], lwd = 1.6)
      }
      points(x, y, col = model_colors[m], pch = model_pch[m],
             bg = model_colors[m], cex = 1.3)
      ## Error bars (IQR)
      arrows(x, lo, x, hi, angle = 90, code = 3,
             length = 0.03, col = model_colors[m], lwd = 0.9)
    }
  }
  
  ## Shared legend in the bottom strip
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = model_labels, col = model_colors[model_order],
         pch = model_pch[model_order], pt.bg = model_colors[model_order],
         lty = model_lty[model_order], lwd = 1.6, ncol = 6, bty = "n",
         cex = 1.0, pt.cex = 1.2)
}


## ============================================================
## Figure 1: Harrell's C-index vs censoring rate
## ============================================================
dev.new(width = 10, height = 4.2)
plot_metric_panels(
  metric_med = "harrell_c_med",
  metric_q25 = "harrell_c_q25",
  metric_q75 = "harrell_c_q75",
  ylab       = "Harrell's C-index",
  main_prefix = "C-index"
)
save_and_show("fig1_c_index.png", width = 10, height = 4.2)


## ============================================================
## Figure 2: Integrated Brier score vs censoring rate
## ============================================================
dev.new(width = 10, height = 4.2)
plot_metric_panels(
  metric_med = "ibs_med",
  metric_q25 = "ibs_q25",
  metric_q75 = "ibs_q75",
  ylab       = "Integrated Brier score",
  main_prefix = "IBS",
  higher_better = FALSE
)
save_and_show("fig2_ibs.png", width = 10, height = 4.2)


## ============================================================
## Figure 3: Calibration slope (grouped bar chart)
## ============================================================
dev.new(width = 10, height = 4.2)

layout_mat <- matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE)
layout(layout_mat, heights = c(5, 1))
par(mar = c(4, 4.2, 2.5, 0.8), family = "serif", cex = 0.85)

for (ds in ds_order) {
  sub <- agg[agg$dataset == ds, ]
  c_vals <- sort(unique(sub$c_obs))
  n_c <- length(c_vals)
  n_m <- length(model_order)
  
  bar_width <- 0.7 / n_m
  group_positions <- seq_len(n_c)
  
  ## Y limits
  ymax <- max(sub$calib_med, na.rm = TRUE) * 1.15
  ymax <- max(ymax, 1.3)
  
  plot(NA, xlim = c(0.3, n_c + 0.7), ylim = c(0, ymax),
       xlab = expression("Censoring rate " * italic(c)),
       ylab = if (ds == ds_order[1]) "Calibration slope" else "",
       main = paste0("(", letters[which(ds_order == ds)], ") ", ds),
       xaxt = "n", las = 1)
  axis(1, at = group_positions, labels = sprintf("%.2f", c_vals))
  grid(col = "#e0e0e0", lty = 1, lwd = 0.5)
  
  ## Reference line at slope = 1
  abline(h = 1.0, lty = 3, col = "black", lwd = 1.0)
  
  for (i in seq_along(model_order)) {
    m <- model_order[i]
    ## Bar position: centered within each group
    offset <- (i - (n_m + 1) / 2) * bar_width
    
    for (g in seq_len(n_c)) {
      row <- sub[sub$model == m & abs(sub$c_obs - c_vals[g]) < 0.02, ]
      if (nrow(row) == 0) next
      val <- row$calib_med[1]
      xpos <- group_positions[g] + offset
      
      rect(xpos - bar_width / 2, 0, xpos + bar_width / 2, val,
           col = model_colors[m], border = "white", lwd = 0.3)
    }
  }
}

## Shared legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = model_labels, fill = model_colors[model_order],
       border = "white", ncol = 6, bty = "n", cex = 1.0)

save_and_show("fig3_calib_slope.png", width = 10, height = 4.2)


## ============================================================
## Figure 4: Empirical partition of the (c, EPV) plane
##
## Each point = one (dataset, regime) cell, labeled by the
## model that achieved the highest median C-index. This is the
## key figure for Section VI.
## ============================================================

## Compute best model per (dataset, regime) cell
best_per_cell <- do.call(rbind, lapply(
  split(agg, interaction(agg$dataset, agg$regime, drop = TRUE)),
  function(g) {
    ## Cell-level c and EPV (median across models)
    c_cell   <- median(g$c_obs, na.rm = TRUE)
    epv_cell <- median(g$epv_obs, na.rm = TRUE)
    
    ## Best model by C-index
    idx_best <- which.max(g$harrell_c_med)
    best     <- g$model[idx_best]
    best_c   <- g$harrell_c_med[idx_best]
    
    ## Second best
    g2 <- g[-idx_best, ]
    idx2 <- which.max(g2$harrell_c_med)
    gap  <- best_c - g2$harrell_c_med[idx2]
    
    data.frame(
      dataset  = g$dataset[1],
      regime   = g$regime[1],
      c_obs    = c_cell,
      epv_obs  = epv_cell,
      best     = best,
      best_c   = best_c,
      gap      = gap,
      stringsAsFactors = FALSE
    )
  }
))
rownames(best_per_cell) <- NULL

dev.new(width = 6.5, height = 5)
par(mar = c(4.5, 4.5, 2, 1), family = "serif", cex = 1.0)

## Point shape and color by winning model
pt_col <- model_colors[best_per_cell$best]
pt_pch <- model_pch[best_per_cell$best]

## Dataset shape for the label
ds_shapes <- c("SUPPORT" = "S", "METABRIC" = "M", "TCGA-BRCA" = "T")

plot(best_per_cell$c_obs, best_per_cell$epv_obs,
     xlim = c(0.25, 0.95), ylim = c(0, max(best_per_cell$epv_obs) * 1.25),
     xlab = expression("Censoring rate " * italic(c)),
     ylab = expression("Events per variable (EPV)"),
     main = "Empirical partition of the (c, EPV) plane",
     pch = pt_pch, bg = pt_col, col = pt_col,
     cex = 2.5, lwd = 1.5, las = 1)
grid(col = "#e0e0e0", lty = 1, lwd = 0.5)

## Re-plot points on top of grid
points(best_per_cell$c_obs, best_per_cell$epv_obs,
       pch = pt_pch, bg = pt_col, col = pt_col,
       cex = 2.5, lwd = 1.5)

## Add dataset + regime labels next to each point
for (i in seq_len(nrow(best_per_cell))) {
  lbl <- paste0(ds_shapes[best_per_cell$dataset[i]], " (",
                best_per_cell$regime[i], ")")
  ## Offset labels to avoid overlap
  xoff <- 0.02
  yoff <- best_per_cell$epv_obs[i] * 0.06
  ## Adjust direction based on position
  if (best_per_cell$c_obs[i] > 0.8) xoff <- -0.07
  text(best_per_cell$c_obs[i] + xoff,
       best_per_cell$epv_obs[i] + yoff,
       labels = lbl, cex = 0.7, adj = 0)
}

## Add annotations showing best model name next to points
for (i in seq_len(nrow(best_per_cell))) {
  m_label <- model_labels[match(best_per_cell$best[i], model_order)]
  yoff2 <- -best_per_cell$epv_obs[i] * 0.08
  xoff2 <- 0.02
  if (best_per_cell$c_obs[i] > 0.8) xoff2 <- -0.07
  text(best_per_cell$c_obs[i] + xoff2,
       best_per_cell$epv_obs[i] + yoff2,
       labels = m_label, cex = 0.65, adj = 0, font = 2,
       col = model_colors[best_per_cell$best[i]])
}

## Legend for model colors
legend("topright",
       legend = model_labels[model_order %in% unique(best_per_cell$best)],
       col = model_colors[model_order[model_order %in% unique(best_per_cell$best)]],
       pch = model_pch[model_order[model_order %in% unique(best_per_cell$best)]],
       pt.bg = model_colors[model_order[model_order %in% unique(best_per_cell$best)]],
       pt.cex = 1.5, bty = "n", cex = 0.85,
       title = "Best model (C-index)")

## Reference line at EPV = 10 (classical rule)
abline(h = 10, lty = 2, col = "#999999", lwd = 1.2)
text(0.27, 10.4, "Classical EPV = 10 guideline", cex = 0.65,
     col = "#999999", adj = 0)

save_and_show("fig4_partition.png", width = 6.5, height = 5)


## ============================================================
## Figure 5: Training time comparison (log scale bar chart)
## ============================================================
dev.new(width = 6, height = 4)
par(mar = c(4, 7, 2, 1), family = "serif", cex = 0.9)

## Pooled median fit time per model
time_pooled <- tapply(res$fit_time_sec, res$model, median, na.rm = TRUE)
time_pooled <- time_pooled[model_order]

barpos <- barplot(log10(time_pooled), horiz = TRUE,
                  names.arg = model_labels, las = 1,
                  col = model_colors[model_order],
                  border = "white",
                  xlab = expression("Median training time per fold (log"[10] * " seconds)"),
                  main = "Computational cost",
                  xlim = c(-2, 3))

## Add actual times as text
for (i in seq_along(time_pooled)) {
  t_val <- time_pooled[i]
  if (t_val >= 1) {
    lbl <- sprintf("%.0f s", t_val)
  } else {
    lbl <- sprintf("%.2f s", t_val)
  }
  text(log10(t_val) + 0.15, barpos[i], lbl, cex = 0.75, adj = 0)
}

save_and_show("fig5_fit_time.png", width = 6, height = 4)


## ============================================================
## Summary: print the partition table for Section VI
## ============================================================
cat("\n=== Partition table (for Table II in the paper) ===\n")
print(best_per_cell[order(best_per_cell$epv_obs), 
                    c("dataset", "regime", "c_obs", "epv_obs",
                      "best", "best_c", "gap")],
      row.names = FALSE)

cat("\nAll figures saved to: ", fig_dir, "\n")
cat("Files:\n")
cat(paste(" ", list.files(fig_dir, pattern = "\\.png$"), collapse = "\n"), "\n")