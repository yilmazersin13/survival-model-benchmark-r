## ============================================================
## 04_cv_scaffolding.R
##
## Cross-validation index generators for the resampling scheme
## of Section IV-D:
##
##   - Outer: 5 x 2 repeated cross-validation. Five independent
##     repetitions of two-fold cross-validation with different
##     random seeds, yielding 10 outer train/test pairs per
##     model-dataset-regime combination.
##
##   - Inner: 5-fold cross-validation on the training portion of
##     each outer fold, used for hyperparameter selection. No
##     information from the outer test fold leaks in.
##
## This module produces ONLY index sets. It does not touch the
## data values, it does not fit any model, it does not compute
## any metric. The fitters in 05_model_fitters.R and the
## evaluation pipeline in main.R consume the index sets and do
## the actual work.
##
## Two structural choices worth flagging:
##
##   1. The outer scheme is stratified by event status. With
##      five repetitions of two-fold CV on small-event datasets
##      such as TCGA-BRCA (n ~ 1000, baseline censoring ~ 0.86,
##      so only ~150 events), unstratified random splits can
##      produce folds with wildly different event counts, which
##      destabilizes the C-index and IBS comparisons. Stratifying
##      on `event` keeps the event proportion approximately
##      constant across folds while preserving the randomization
##      that justifies 5x2 CV. This is standard practice in
##      survival benchmarks; Section IV-D does not forbid it and
##      the description in the paper is consistent with it.
##
##   2. Each outer repetition uses its own seed, derived from a
##      single base seed. This means the entire 5x2 structure
##      for a given dataset-regime combination is reproducible
##      from one integer, which simplifies the main script.
## ============================================================

## ------------------------------------------------------------
## Internal helper: stratified split of indices into k folds.
##
## Given a binary stratification vector `strat` (typically the
## event indicator) and a number of folds k, returns an integer
## vector of fold assignments of length n, with each stratum
## distributed approximately evenly across folds.
## ------------------------------------------------------------
.stratified_fold_assign <- function(strat, k, seed) {
  n <- length(strat)
  fold_id <- integer(n)
  
  set.seed(seed)
  for (lvl in unique(strat)) {
    idx <- which(strat == lvl)
    n_lvl <- length(idx)
    ## Distribute indices in this stratum across k folds as evenly
    ## as possible, then shuffle so that the assignment is random
    ## within the stratum.
    assign_lvl <- rep(seq_len(k), length.out = n_lvl)
    assign_lvl <- sample(assign_lvl)
    fold_id[idx] <- assign_lvl
  }
  fold_id
}


## ============================================================
## make_outer_cv()
##
## Generate the 5 x 2 outer cross-validation index structure
## for a given dataset, stratified by event status.
##
## Arguments:
##   event       integer vector in {0,1}, length n. The event
##               indicator delta_i (or delta_i^(c) if augmentation
##               has already been applied to the dataset).
##   n_reps      integer, number of repetitions; default 5
##   n_folds     integer, number of folds per repetition; default 2
##   base_seed   integer, master seed; per-repetition seeds are
##               derived deterministically from this
##
## Returns a list of length n_reps * n_folds (default 10), each
## element of the form:
##   list(rep = r, fold = k, train = <indices>, test = <indices>)
##
## The train and test indices within a single (rep, fold) entry
## are disjoint and together cover all rows. Across the two folds
## of one repetition, the train/test roles are exchanged.
## ============================================================
make_outer_cv <- function(event, n_reps = 5, n_folds = 2,
                          base_seed = 20260413) {
  if (!all(event %in% c(0L, 1L, 0, 1))) {
    stop("event must be 0/1")
  }
  n <- length(event)
  all_idx <- seq_len(n)
  
  splits <- vector("list", n_reps * n_folds)
  pos <- 1L
  
  for (r in seq_len(n_reps)) {
    ## Derive a per-repetition seed from the base seed. Using
    ## base_seed + r * <large prime> gives well-separated streams
    ## without requiring the L'Ecuyer machinery, which is overkill
    ## here because we are not running things in parallel inside
    ## the CV loop.
    rep_seed <- base_seed + r * 100003L
    
    fold_id <- .stratified_fold_assign(event, k = n_folds, seed = rep_seed)
    
    for (k in seq_len(n_folds)) {
      test_idx  <- all_idx[fold_id == k]
      train_idx <- all_idx[fold_id != k]
      
      splits[[pos]] <- list(
        rep   = r,
        fold  = k,
        train = train_idx,
        test  = test_idx
      )
      pos <- pos + 1L
    }
  }
  
  ## Sanity check: every row should appear in `test` exactly
  ## n_reps times across the full structure.
  test_counts <- tabulate(unlist(lapply(splits, `[[`, "test")), nbins = n)
  if (any(test_counts != n_reps)) {
    stop("make_outer_cv: internal error, test coverage is inconsistent")
  }
  
  splits
}


## ============================================================
## make_inner_cv()
##
## Generate the inner k-fold cross-validation index structure
## for hyperparameter tuning within a single outer training
## fold. Stratified by event status.
##
## Arguments:
##   train_idx   integer vector, the row indices of the OUTER
##               training fold (these are positions in the
##               original n-row dataset, NOT 1..length(train_idx))
##   event_full  integer vector in {0,1}, length n, the event
##               indicator of the full dataset. The function
##               extracts event_full[train_idx] for stratification.
##   n_folds     integer, number of inner folds; default 5
##   seed        integer, RNG seed for the inner split
##
## Returns a list of length n_folds, each element of the form:
##   list(fold = k, train = <indices>, valid = <indices>)
##
## The returned indices are positions in the ORIGINAL dataset,
## not positions within `train_idx`. This is the convention used
## throughout the rest of the library: every index is always a
## row of the original n-row dataset, regardless of how deeply
## nested the resampling is. It avoids the very common bug of
## indexing a subset by a sub-subset's local positions.
## ============================================================
make_inner_cv <- function(train_idx, event_full, n_folds = 5,
                          seed = 1L) {
  if (length(train_idx) < n_folds) {
    stop(sprintf("make_inner_cv: training fold has only %d rows, ",
                 length(train_idx)),
         sprintf("cannot make %d inner folds", n_folds))
  }
  if (!all(event_full[train_idx] %in% c(0L, 1L, 0, 1))) {
    stop("event_full must be 0/1")
  }
  
  strat_local <- event_full[train_idx]
  fold_id <- .stratified_fold_assign(strat_local, k = n_folds, seed = seed)
  
  splits <- vector("list", n_folds)
  for (k in seq_len(n_folds)) {
    valid_local <- which(fold_id == k)
    train_local <- which(fold_id != k)
    
    ## Translate local positions back to original dataset row indices.
    splits[[k]] <- list(
      fold  = k,
      train = train_idx[train_local],
      valid = train_idx[valid_local]
    )
  }
  splits
}


## ============================================================
## summarize_cv_structure()
##
## Diagnostic helper. Given an outer CV object and the full event
## vector, print a small table showing how many training events,
## training non-events, test events, and test non-events fall in
## each (rep, fold) cell. Used at the start of the main script to
## confirm the splits look sensible before any model is fit; not
## required for the benchmark itself.
## ============================================================
summarize_cv_structure <- function(outer_cv, event_full) {
  rows <- lapply(outer_cv, function(s) {
    data.frame(
      rep            = s$rep,
      fold           = s$fold,
      n_train        = length(s$train),
      n_train_events = sum(event_full[s$train] == 1),
      n_test         = length(s$test),
      n_test_events  = sum(event_full[s$test] == 1)
    )
  })
  df <- do.call(rbind, rows)
  df$pct_train_events <- round(100 * df$n_train_events / df$n_train, 1)
  df$pct_test_events  <- round(100 * df$n_test_events  / df$n_test,  1)
  df
}