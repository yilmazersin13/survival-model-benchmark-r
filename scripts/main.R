## ============================================================
## main.R
##
## Orchestration script for the benchmark of Section IV. Loops
## over (dataset x regime x outer fold x model), calls each
## fitter through the uniform interface from
## 05_model_fitters.R / 05b_pycox_wrappers.R, computes metrics
## via 06_evaluation_metrics.R and diagnostics via
## 07_diagnostics.R, and incrementally appends rows to the long
## results table that 08_aggregation_and_stats.R consumes.
##
## Design priorities, in order:
##
##   1. Crash resilience. A 6-10 hour overnight run must not
##      lose all work if one fitter crashes on one fold. Every
##      (cell, model) attempt is wrapped in tryCatch and the
##      results table is appended to disk after every fold.
##      Restarting the script picks up where it left off by
##      reading the existing results CSV.
##
##   2. Honest accounting of failures. When a fit fails or
##      returns NA metrics, that fact is recorded in the row
##      rather than silently dropped. Section V can then report
##      "Cox failed to converge in 10/10 folds on METABRIC" as
##      a real finding.
##
##   3. Fold-aware everything. Preprocessing, augmentation, and
##      EPV are all computed on the training fold of each outer
##      split, never on the full dataset and never on test
##      rows.
##
## Usage:
##   source("main.R")
##   run_benchmark(out_dir = "benchmark_results")
##
## To resume an interrupted run:
##   run_benchmark(out_dir = "benchmark_results", resume = TRUE)
##
## To restrict the run for testing (e.g. one dataset, one
## regime, two folds):
##   run_benchmark(
##     datasets_to_run = "SUPPORT",
##     regimes_to_run  = "baseline",
##     n_outer_folds   = 2,
##     n_random_search = 2
##   )
## ============================================================

## ------------------------------------------------------------
## Source all eight modules. Order matters slightly (modules
## with no dependencies first), but most are independent so
## this is mostly bookkeeping.
## ------------------------------------------------------------
source("01_data_loaders.R")
source("02_preprocessing.R")
source("03_censoring_augmentation.R")
source("04_cv_scaffolding.R")
source("05_model_fitters.R")
source("05b_pycox_wrappers.R")
source("06_evaluation_metrics.R")
source("07_diagnostics.R")
source("08_aggregation_and_stats.R")


## ============================================================
## Configuration
## ============================================================

## The six models in the benchmark, in the order they will be
## fit per cell. Order is purely cosmetic except that BEN-Cox is
## last because it is the slowest and we want the cheaper models
## to fail fast if something is wrong with the cell.
.MODELS <- c("cox", "lasso_cox", "rsf", "deepsurv", "cox_time", "ben_cox")

## Map model names to their fitter functions. Used by the per-
## cell fit loop so the dispatch is data, not a switch.
.FITTERS <- list(
  cox       = fit_cox,
  lasso_cox = fit_lasso_cox,
  rsf       = fit_rsf,
  deepsurv  = fit_deepsurv,
  cox_time  = fit_cox_time,
  ben_cox   = fit_ben_cox
)

## Censoring regimes in the order they will be evaluated.
## "baseline" means use the dataset's native censoring. The
## numeric values are the targets passed to augment_dataset();
## values at or below the dataset's baseline are no-ops and
## the row is recorded only once (under "baseline").
.REGIMES <- c("baseline", "0.30", "0.50", "0.70")


## ------------------------------------------------------------
## Logging helper: print a timestamped message and flush.
## ------------------------------------------------------------
.log <- function(...) {
  msg <- sprintf("[%s] %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n", sep = "")
  flush.console()
}


## ------------------------------------------------------------
## Build the per-fold results row in the canonical schema that
## 08_aggregation_and_stats.R expects. Centralized here so the
## column structure stays consistent across success and failure
## paths.
## ------------------------------------------------------------
.make_row <- function(dataset, regime, model_name, fold, rep,
                      c_target, c_achieved,
                      n_events,
                      harrell_c = NA_real_, uno_c = NA_real_,
                      ibs = NA_real_, calib_slope = NA_real_,
                      epv = NA_real_, epv_raw = NA_real_,
                      p_eff = NA_integer_, p_raw = NA_integer_,
                      fit_time_sec = NA_real_,
                      converged = NA, fit_status = "ok",
                      ph_global_p = NA_real_) {
  data.frame(
    dataset      = dataset,
    regime       = regime,
    model        = model_name,
    fold         = as.integer(fold),
    rep          = as.integer(rep),
    c_target     = c_target,
    c_achieved   = c_achieved,
    n_events     = as.integer(n_events),
    harrell_c    = harrell_c,
    uno_c        = uno_c,
    ibs          = ibs,
    calib_slope  = calib_slope,
    epv          = epv,
    epv_raw      = epv_raw,
    p_eff        = as.integer(p_eff),
    p_raw        = as.integer(p_raw),
    fit_time_sec = fit_time_sec,
    converged    = converged,
    fit_status   = fit_status,
    ph_global_p  = ph_global_p,
    stringsAsFactors = FALSE
  )
}


## ------------------------------------------------------------
## Append one row to the results CSV. If the file does not
## exist, write with header; otherwise append without header.
## Done one row at a time so the run is crash-resilient: even
## if R is killed mid-fold, the rows already written are safe.
## ------------------------------------------------------------
.append_row <- function(row, csv_path) {
  if (file.exists(csv_path)) {
    utils::write.table(row, csv_path, sep = ",",
                       row.names = FALSE, col.names = FALSE,
                       append = TRUE, qmethod = "double")
  } else {
    utils::write.csv(row, csv_path, row.names = FALSE)
  }
}


## ------------------------------------------------------------
## Read existing results CSV (for resume mode) and return a
## set of (dataset, regime, fold, model) keys already complete.
## ------------------------------------------------------------
.completed_keys <- function(csv_path) {
  if (!file.exists(csv_path)) return(character(0))
  done <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (is.null(done) || nrow(done) == 0) return(character(0))
  paste(done$dataset, done$regime, done$fold, done$model, sep = "@@")
}


## ============================================================
## run_one_cell()
##
## Fit all six models on one (dataset, regime, fold) cell and
## append rows to the results CSV. This is the inner unit of
## work; everything else is loops around it.
##
## Arguments:
##   dataset_obj      a dataset object from 01_data_loaders.R
##                    (already augmented to the regime, if any)
##   regime           character label for the regime
##   c_target         numeric target c (NA for baseline)
##   c_achieved       numeric achieved c after augmentation
##   outer_split      one element of the list from
##                    make_outer_cv()
##   csv_path         path to the results CSV
##   already_done     character vector of completed keys
##                    (skip work where key matches)
##   horizon_quantile numeric in (0,1); the per-dataset
##                    horizon for calibration is computed as
##                    this quantile of training event times
##   eval_quantiles   numeric vector in (0,1); IBS evaluation
##                    grid as quantiles of training event times
##   models_to_run    character vector subset of .MODELS
##   n_random_search  passed to deepsurv/cox_time fitters
## ============================================================
run_one_cell <- function(dataset_obj, regime, c_target, c_achieved,
                         outer_split, csv_path,
                         already_done = character(0),
                         horizon_quantile = 0.75,
                         eval_quantiles   = seq(0.1, 0.85, by = 0.1),
                         models_to_run    = .MODELS,
                         n_random_search  = 8,
                         seed_offset      = 0L) {
  ds_name <- dataset_obj$dataset_name
  fold    <- outer_split$fold
  rep     <- outer_split$rep
  cell_id <- sprintf("%s @ %s @ rep%d/fold%d",
                     ds_name, regime, rep, fold)
  
  ## --- Preprocess the fold using the (possibly augmented) dataset ---
  pp <- tryCatch(
    preprocess_fold(dataset_obj,
                    train_idx = outer_split$train,
                    test_idx  = outer_split$test),
    error = function(e) {
      .log("  ! preprocessing failed for ", cell_id, ": ",
           conditionMessage(e))
      NULL
    }
  )
  if (is.null(pp)) return(invisible(NULL))
  
  p_raw    <- ncol(pp$X_train)
  n_events <- sum(pp$event_train == 1)
  
  ## --- Per-cell PH diagnostic on the unpenalized Cox ---
  ## Run once per cell, then attach the same global p-value to
  ## every model's row. The diagnostic is a property of the
  ## DATASET cell, not of any individual fitter.
  ph <- tryCatch(
    ph_test_global(pp$X_train, pp$time_train, pp$event_train),
    error = function(e) list(global_p = NA_real_, converged = FALSE)
  )
  
  ## --- Per-cell horizon and IBS evaluation grid ---
  event_times_train <- pp$time_train[pp$event_train == 1]
  if (length(event_times_train) < 5) {
    .log("  ! cell ", cell_id, " has fewer than 5 training events; skipping")
    return(invisible(NULL))
  }
  horizon    <- as.numeric(quantile(event_times_train, horizon_quantile))
  eval_times <- as.numeric(quantile(event_times_train, eval_quantiles))
  eval_times <- sort(unique(eval_times))
  
  ## --- Loop over the six fitters ---
  for (m in models_to_run) {
    key <- paste(ds_name, regime, fold, m, sep = "@@")
    if (key %in% already_done) {
      .log("  - skip (already done): ", cell_id, " :: ", m)
      next
    }
    
    fitter <- .FITTERS[[m]]
    if (is.null(fitter)) {
      .log("  ! no fitter for model: ", m)
      next
    }
    
    .log("  > fit ", cell_id, " :: ", m)
    
    ## Per-model seed: combine cell identity with model name
    ## for reproducibility within a fitter, while still varying
    ## across cells.
    cell_seed <- 100000L * fold + 1000L * rep +
      utils::adist(m, "ben_cox")[1, 1] + seed_offset
    
    fit_result <- tryCatch(
      {
        if (m %in% c("deepsurv", "cox_time")) {
          fitter(pp$X_train, pp$time_train, pp$event_train,
                 seed = cell_seed,
                 n_random_search = n_random_search)
        } else {
          fitter(pp$X_train, pp$time_train, pp$event_train,
                 seed = cell_seed)
        }
      },
      error = function(e) {
        .log("    ! fit error: ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(fit_result)) {
      ## Fit hard-failed; record an honest failure row so the
      ## aggregation step knows this cell-model combination was
      ## attempted.
      row <- .make_row(
        dataset = ds_name, regime = regime, model_name = m,
        fold = fold, rep = rep,
        c_target = c_target, c_achieved = c_achieved,
        n_events = n_events, p_raw = p_raw,
        fit_status = "fit_error",
        ph_global_p = ph$global_p
      )
      .append_row(row, csv_path)
      next
    }
    
    ## --- Compute metrics on the test fold ---
    metrics <- tryCatch(
      {
        risk <- fit_result$predict_risk(pp$X_test)
        surv_matrix <- fit_result$predict_survival(pp$X_test, eval_times)
        surv_at_h   <- as.numeric(
          fit_result$predict_survival(pp$X_test, horizon))
        
        list(
          harrell_c   = harrell_c(pp$time_test, pp$event_test, risk),
          uno_c       = uno_c(pp$time_test, pp$event_test, risk,
                              time_train  = pp$time_train,
                              event_train = pp$event_train),
          ibs         = integrated_brier_score(
            pp$time_test, pp$event_test, surv_matrix,
            eval_times  = eval_times,
            time_train  = pp$time_train,
            event_train = pp$event_train),
          calib_slope = calibration_slope(
            pp$time_test, pp$event_test, surv_at_h,
            horizon = horizon)
        )
      },
      error = function(e) {
        .log("    ! metric error: ", conditionMessage(e))
        list(harrell_c = NA_real_, uno_c = NA_real_,
             ibs = NA_real_, calib_slope = NA_real_)
      }
    )
    
    ## --- EPV (model-specific) ---
    epv_info <- tryCatch(
      compute_epv(pp$event_train, fit_result, p_raw = p_raw),
      error = function(e) {
        list(epv = NA_real_, epv_raw = NA_real_,
             p_eff = NA_integer_, n_events = n_events)
      }
    )
    
    converged <- if (!is.null(fit_result$hyperparams$converged)) {
      fit_result$hyperparams$converged
    } else {
      TRUE
    }
    
    row <- .make_row(
      dataset      = ds_name, regime = regime, model_name = m,
      fold         = fold, rep = rep,
      c_target     = c_target, c_achieved = c_achieved,
      n_events     = n_events,
      harrell_c    = metrics$harrell_c,
      uno_c        = metrics$uno_c,
      ibs          = metrics$ibs,
      calib_slope  = metrics$calib_slope,
      epv          = epv_info$epv,
      epv_raw      = epv_info$epv_raw,
      p_eff        = epv_info$p_eff,
      p_raw        = p_raw,
      fit_time_sec = fit_result$fit_time_sec,
      converged    = converged,
      fit_status   = "ok",
      ph_global_p  = ph$global_p
    )
    .append_row(row, csv_path)
    
    ## --- Aggressive memory cleanup ---
    ## BEN-Cox carries a Stan model object that can be hundreds
    ## of MB. DeepSurv/Cox-Time carry PyTorch model objects.
    ## Clearing fit_result and forcing a GC between models keeps
    ## peak memory bounded over the long run.
    rm(fit_result)
    gc(verbose = FALSE, full = TRUE)
  }
  
  invisible(NULL)
}


## ============================================================
## run_benchmark()
##
## Top-level orchestrator. Iterates over datasets x regimes x
## outer folds, builds the augmented dataset and CV indices for
## each cell, and calls run_one_cell() to do the actual work.
##
## Arguments are documented at the top of this file; defaults
## reproduce Section IV verbatim.
## ============================================================
run_benchmark <- function(
    out_dir          = "benchmark_results",
    datasets_to_run  = c("METABRIC", "SUPPORT", "TCGA-BRCA"),
    regimes_to_run   = .REGIMES,
    models_to_run    = .MODELS,
    n_outer_folds    = 10L,
    n_random_search  = 8L,
    base_seed        = 20260413L,
    resume           = FALSE,
    metabric_file    = "METABRIC_RNA_Mutation.csv",
    support_file     = "support.tsv",
    tcga_file        = "tcga_brca_data.tsv",
    tcga_expression_file = NULL
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  csv_path <- file.path(out_dir, "raw_results.csv")
  log_path <- file.path(out_dir, "run.log")
  
  ## Tee log output to a file as well as the console.
  log_con <- file(log_path, open = "at")
  sink(log_con, append = TRUE, type = "output", split = TRUE)
  on.exit({
    sink(type = "output")
    close(log_con)
  })
  
  if (!resume && file.exists(csv_path)) {
    bak <- file.path(out_dir,
                     sprintf("raw_results_backup_%s.csv",
                             format(Sys.time(), "%Y%m%d_%H%M%S")))
    file.rename(csv_path, bak)
    .log("Existing results backed up to ", bak)
  }
  
  already_done <- if (resume) .completed_keys(csv_path) else character(0)
  if (length(already_done) > 0) {
    .log("Resuming run: ", length(already_done), " cells already complete")
  }
  
  .log("=== Benchmark started ===")
  .log("Datasets: ", paste(datasets_to_run, collapse = ", "))
  .log("Regimes:  ", paste(regimes_to_run,  collapse = ", "))
  .log("Models:   ", paste(models_to_run,   collapse = ", "))
  .log("Outer folds: ", n_outer_folds, "; deep-model search size: ",
       n_random_search)
  
  ## --- Load datasets once ---
  all_data <- list()
  if ("METABRIC"   %in% datasets_to_run) all_data$METABRIC   <- load_metabric(metabric_file)
  if ("SUPPORT"    %in% datasets_to_run) all_data$SUPPORT    <- load_support(support_file)
  if ("TCGA-BRCA"  %in% datasets_to_run) all_data$`TCGA-BRCA` <- load_tcga_brca(tcga_file, tcga_expression_file)
  
  t0_overall <- Sys.time()
  
  for (ds_name in datasets_to_run) {
    ds <- all_data[[ds_name]]
    if (is.null(ds)) {
      .log("! dataset not loaded: ", ds_name, "; skipping")
      next
    }
    
    for (regime in regimes_to_run) {
      ## --- Build the augmented dataset for this regime ---
      if (regime == "baseline") {
        ds_aug <- ds
        ds_aug$augmentation <- list(c_target = NA_real_,
                                    c_achieved = ds$baseline_censoring,
                                    augmented = FALSE)
      } else {
        c_target <- as.numeric(regime)
        if (c_target <= ds$baseline_censoring + 0.01) {
          .log("- skip ", ds_name, " @ regime ", regime,
               " (baseline ", round(ds$baseline_censoring, 3),
               " already at or above target)")
          next
        }
        ## Note: augment_dataset uses a single seed; the per-fold
        ## CV indices come from a separate seeded generator, so
        ## the fold structure is reproducible independently of
        ## the augmentation draws.
        ds_aug <- augment_dataset(ds, c_target = c_target,
                                  seed = base_seed + 7L)
      }
      
      c_target   <- ds_aug$augmentation$c_target
      c_achieved <- ds_aug$augmentation$c_achieved
      
      ## --- Build outer CV on the augmented event vector ---
      outer_cv <- make_outer_cv(ds_aug$event,
                                base_seed = base_seed)
      ## Truncate to n_outer_folds for testing runs.
      outer_cv <- outer_cv[seq_len(min(length(outer_cv), n_outer_folds))]
      
      .log("=== ", ds_name, " @ regime ", regime,
           " (c_achieved = ", round(c_achieved, 3),
           ", n_events = ", sum(ds_aug$event == 1),
           ", folds = ", length(outer_cv), ") ===")
      
      for (i in seq_along(outer_cv)) {
        run_one_cell(
          dataset_obj      = ds_aug,
          regime           = regime,
          c_target         = c_target,
          c_achieved       = c_achieved,
          outer_split      = outer_cv[[i]],
          csv_path         = csv_path,
          already_done     = already_done,
          models_to_run    = models_to_run,
          n_random_search  = n_random_search
        )
      }
    }
  }
  
  total_min <- as.numeric(Sys.time() - t0_overall, units = "mins")
  .log("=== Benchmark finished in ", round(total_min, 1), " minutes ===")
  
  ## --- Final aggregation, stats, and partition ---
  .log("Computing aggregations and writing summary tables ...")
  results_df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  paths <- save_results(results_df, out_dir = out_dir)
  .log("Outputs:")
  for (p in paths) .log("  ", p)
  
  invisible(results_df)
}


## ============================================================
## Entry point
## ============================================================
##
## Interactive use: source this file then call run_benchmark().
## Examples:
##
##   # Full benchmark (default settings, ~6-10 hours)
   run_benchmark()
##
##   # Smoke test: SUPPORT, baseline regime only, two folds,
##   # tiny deep-model search. Should complete in 5-10 minutes.
#   run_benchmark(
#     datasets_to_run = "SUPPORT",
#     regimes_to_run  = "baseline",
#     n_outer_folds   = 2,
#     n_random_search = 2,
#     out_dir         = "test_run"
#   )
##
##   # Resume an interrupted run
##   run_benchmark(resume = TRUE)
##
## ============================================================