
##
##   dataset       character (e.g. "METABRIC")
##   regime        character or numeric (e.g. "0.50" or 0.50)
##   model         character (one of the six model_class strings)
##   fold          integer (1..10 across the 5x2 outer CV)
##   rep           integer (1..5)
##   harrell_c     numeric, possibly NA
##   uno_c         numeric, possibly NA
##   ibs           numeric, possibly NA
##   calib_slope   numeric, possibly NA
##   epv           numeric, model-specific EPV from compute_epv()
##   epv_raw       numeric, raw events-to-p ratio
##   p_eff         integer
##   n_events      integer
##   fit_time_sec  numeric
##   converged     logical or NA
##

.summarize_vec <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(c(median = NA_real_, q25 = NA_real_, q75 = NA_real_, n = 0L))
  }
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE,
                 na.rm = TRUE, type = 7)
  c(median = qs[2], q25 = qs[1], q75 = qs[3], n = as.integer(n))
}


aggregate_metrics <- function(results_df,
                              metrics = c("harrell_c", "uno_c", "ibs",
                                          "calib_slope",
                                          "epv", "epv_raw", "p_eff",
                                          "n_events", "fit_time_sec")) {
  required <- c("dataset", "regime", "model", "fold")
  missing <- setdiff(required, names(results_df))
  if (length(missing) > 0) {
    stop(sprintf("aggregate_metrics: results_df missing columns: %s",
                 paste(missing, collapse = ", ")))
  }
  
  ## Group by (dataset, regime, model). Use split() over an
  ## interaction factor rather than aggregate() because we want
  ## to compute multiple statistics per group with custom names,
  ## which aggregate() does not support cleanly.
  group_key <- interaction(results_df$dataset,
                           results_df$regime,
                           results_df$model,
                           drop = TRUE, sep = "@@")
  groups <- split(results_df, group_key)
  
  rows <- lapply(groups, function(g) {
    out <- data.frame(
      dataset = g$dataset[1],
      regime  = g$regime[1],
      model   = g$model[1],
      stringsAsFactors = FALSE
    )
    for (m in metrics) {
      if (!m %in% names(g)) next
      stats <- .summarize_vec(g[[m]])
      out[[paste0(m, "_median")]] <- stats[["median"]]
      out[[paste0(m, "_q25")]]    <- stats[["q25"]]
      out[[paste0(m, "_q75")]]    <- stats[["q75"]]
      out[[paste0(m, "_n")]]      <- stats[["n"]]
    }
    out
  })
  
  do.call(rbind, rows)
}



pairwise_wilcoxon_holm <- function(results_df,
                                   metrics = c("harrell_c", "uno_c",
                                               "ibs", "calib_slope"),
                                   alternative = "two.sided") {
  required <- c("dataset", "regime", "model", "fold")
  missing <- setdiff(required, names(results_df))
  if (length(missing) > 0) {
    stop(sprintf("pairwise_wilcoxon_holm: missing columns: %s",
                 paste(missing, collapse = ", ")))
  }
  
  rows <- list()
  row_idx <- 1L
  
  ## Iterate over (dataset, regime) cells.
  cell_key <- interaction(results_df$dataset, results_df$regime,
                          drop = TRUE, sep = "@@")
  cells <- split(results_df, cell_key)
  
  for (cell in cells) {
    ds  <- cell$dataset[1]
    reg <- cell$regime[1]
    models_here <- sort(unique(cell$model))
    if (length(models_here) < 2) next
    
    ## All unordered pairs of models within this cell.
    pairs <- utils::combn(models_here, 2, simplify = FALSE)
    
    for (metric in metrics) {
      if (!metric %in% names(cell)) next
      
      ## Collect raw p-values for Holm correction within this
      ## (dataset, regime, metric) stratum. We accumulate first,
      ## then adjust, then write rows.
      raw <- lapply(pairs, function(pr) {
        a <- pr[1]; b <- pr[2]
        ga <- cell[cell$model == a, ]
        gb <- cell[cell$model == b, ]
        ## Pair on fold so test is paired, not pooled.
        merged <- merge(ga[, c("fold", metric), drop = FALSE],
                        gb[, c("fold", metric), drop = FALSE],
                        by = "fold", suffixes = c("_a", "_b"))
        va <- merged[[paste0(metric, "_a")]]
        vb <- merged[[paste0(metric, "_b")]]
        keep <- is.finite(va) & is.finite(vb)
        va <- va[keep]; vb <- vb[keep]
        n_paired <- length(va)
        
        if (n_paired < 5) {
          return(list(
            model_a = a, model_b = b,
            median_a = if (length(va) > 0) median(va) else NA_real_,
            median_b = if (length(vb) > 0) median(vb) else NA_real_,
            median_diff = NA_real_,
            wilcoxon_p = NA_real_,
            n_paired = n_paired
          ))
        }
        
        wt <- tryCatch(
          suppressWarnings(
            wilcox.test(va, vb, paired = TRUE,
                        alternative = alternative,
                        exact = FALSE)
          ),
          error = function(e) NULL
        )
        list(
          model_a = a, model_b = b,
          median_a = median(va),
          median_b = median(vb),
          median_diff = median(va - vb),
          wilcoxon_p = if (is.null(wt)) NA_real_ else unname(wt$p.value),
          n_paired = n_paired
        )
      })
      
      raw_p <- vapply(raw, function(r) r$wilcoxon_p, numeric(1))
      ## Holm correction within this stratum. NA p-values stay NA.
      p_adj <- rep(NA_real_, length(raw_p))
      ok <- is.finite(raw_p)
      if (any(ok)) p_adj[ok] <- p.adjust(raw_p[ok], method = "holm")
      
      for (k in seq_along(raw)) {
        r <- raw[[k]]
        conclusion <- if (!is.finite(p_adj[k])) {
          "insufficient_data"
        } else if (p_adj[k] >= 0.05) {
          "ns"
        } else if (is.finite(r$median_diff) && r$median_diff > 0) {
          "a > b"
        } else if (is.finite(r$median_diff) && r$median_diff < 0) {
          "b > a"
        } else {
          "ns"
        }
        
        rows[[row_idx]] <- data.frame(
          dataset       = ds,
          regime        = reg,
          metric        = metric,
          model_a       = r$model_a,
          model_b       = r$model_b,
          median_a      = r$median_a,
          median_b      = r$median_b,
          median_diff   = r$median_diff,
          wilcoxon_p    = r$wilcoxon_p,
          p_holm        = p_adj[k],
          n_paired      = r$n_paired,
          conclusion    = conclusion,
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }
  
  if (length(rows) == 0) {
    return(data.frame(
      dataset = character(0), regime = character(0),
      metric = character(0),
      model_a = character(0), model_b = character(0),
      median_a = numeric(0), median_b = numeric(0),
      median_diff = numeric(0),
      wilcoxon_p = numeric(0), p_holm = numeric(0),
      n_paired = integer(0), conclusion = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

construct_partition <- function(results_df,
                                metric = "harrell_c",
                                higher_is_better = TRUE,
                                epv_column = "epv") {
  required <- c("dataset", "regime", "model", metric, epv_column)
  missing <- setdiff(required, names(results_df))
  if (length(missing) > 0) {
    stop(sprintf("construct_partition: missing columns: %s",
                 paste(missing, collapse = ", ")))
  }
  
  cell_key <- interaction(results_df$dataset, results_df$regime,
                          drop = TRUE, sep = "@@")
  cells <- split(results_df, cell_key)
  
  rows <- lapply(cells, function(cell) {
    ds  <- cell$dataset[1]
    reg <- cell$regime[1]
    
    ## Cell-level c (censoring rate). The benchmark assigns each
    ## (dataset, regime) cell a single target c, but the achieved
    ## c per fold can differ slightly (the bisection tolerance is
    ## 0.01). Report the median achieved c for accuracy.
    c_obs <- if ("c_achieved" %in% names(cell)) {
      median(cell$c_achieved, na.rm = TRUE)
    } else {
      ## Fall back to using the regime label as a numeric.
      suppressWarnings(as.numeric(reg))
    }
    
    ## EPV is reported as the median across folds AND across
    ## models. For cells where p_eff differs by model (lasso vs
    ## ben_cox vs cox), this produces a cell-level EPV summary
    ## that is interpretable on the partition plot.
    epv_obs <- median(cell[[epv_column]], na.rm = TRUE)
    
    ## Per-model median of the chosen metric.
    by_model <- split(cell[[metric]], cell$model)
    medians <- vapply(by_model, function(v) {
      vv <- v[is.finite(v)]
      if (length(vv) == 0) NA_real_ else median(vv)
    }, numeric(1))
    medians <- medians[is.finite(medians)]
    
    if (length(medians) == 0) {
      return(data.frame(
        dataset = ds, regime = reg,
        c_observed = c_obs, epv_observed = epv_obs,
        best_model = NA_character_,
        best_metric_median = NA_real_,
        second_model = NA_character_,
        second_metric_median = NA_real_,
        gap_to_second = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    ord <- order(medians, decreasing = higher_is_better)
    ranked <- medians[ord]
    
    best_model <- names(ranked)[1]
    best_med   <- unname(ranked[1])
    second_model <- if (length(ranked) >= 2) names(ranked)[2] else NA_character_
    second_med   <- if (length(ranked) >= 2) unname(ranked[2]) else NA_real_
    gap <- if (is.finite(second_med)) abs(best_med - second_med) else NA_real_
    
    data.frame(
      dataset = ds, regime = reg,
      c_observed = c_obs, epv_observed = epv_obs,
      best_model = best_model,
      best_metric_median = best_med,
      second_model = second_model,
      second_metric_median = second_med,
      gap_to_second = gap,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}

save_results <- function(results_df, out_dir = "benchmark_results") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  paths <- list()
  
  paths$raw <- file.path(out_dir, "raw_results.csv")
  utils::write.csv(results_df, paths$raw, row.names = FALSE)
  
  agg <- aggregate_metrics(results_df)
  paths$agg <- file.path(out_dir, "aggregated.csv")
  utils::write.csv(agg, paths$agg, row.names = FALSE)
  
  pw <- pairwise_wilcoxon_holm(results_df)
  paths$pairwise <- file.path(out_dir, "pairwise_wilcoxon.csv")
  utils::write.csv(pw, paths$pairwise, row.names = FALSE)
  
  part_eff <- construct_partition(results_df, epv_column = "epv")
  paths$part_eff <- file.path(out_dir, "partition.csv")
  utils::write.csv(part_eff, paths$part_eff, row.names = FALSE)
  
  if ("epv_raw" %in% names(results_df)) {
    part_raw <- construct_partition(results_df, epv_column = "epv_raw")
    paths$part_raw <- file.path(out_dir, "partition_raw_epv.csv")
    utils::write.csv(part_raw, paths$part_raw, row.names = FALSE)
  }
  
  message(sprintf("Wrote %d result files to %s", length(paths), out_dir))
  invisible(paths)
}
