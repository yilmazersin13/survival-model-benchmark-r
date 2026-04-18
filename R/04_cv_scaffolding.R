
## 04_cv_scaffolding.R

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
