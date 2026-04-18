## ============================================================
## 05b_pycox_wrappers.R
##
## Wrappers around pycox's deep survival models, accessed from
## R via reticulate. Two fitters:
##
##   fit_deepsurv()   DeepSurv (Katzman et al. 2018), implemented
##                    as pycox.models.CoxPH. Section III-B.
##   fit_cox_time()   Cox-Time (Kvamme, Borgan & Scheel 2019),
##                    implemented as pycox.models.CoxTime.
##                    Section III-B.
##
## Both fitters expose the same eight-slot return list as the
## R-native fitters in 05_model_fitters.R, so downstream code
## can call them interchangeably:
##
##   $model_class       "deepsurv" | "cox_time"
##   $raw_fit           the underlying pycox model object (Python
##                      reference, not directly usable from R)
##   $coef              NULL (deep models do not have linear coef)
##   $predict_risk      function(X_new) -> numeric vector
##   $predict_survival  function(X_new, times) -> matrix
##   $hyperparams       list of final hyperparameters used
##   $fit_time_sec      numeric, wall-clock seconds for the fit
##
## Hyperparameter search. Section IV-D specifies a small random
## search over learning rate, dropout, hidden layer width, and
## number of hidden layers, with early stopping on an inner
## validation split. The wrappers below implement this with a
## fixed budget controlled by the n_random_search argument.
## Defaults are kept modest to match the budget-matched design
## of Section IV-D.
##
## CPU vs GPU. The wrappers are CPU-only by default. PyTorch will
## use whatever device is available; on a machine without CUDA
## this is automatic. No code changes needed for GPU machines.
##
## Reproducibility. The fitters seed both R's RNG (for the
## random search and the validation split) and Python's torch
## RNG (for network weight initialization and the training
## loop). Given the same seed, the same data, and the same
## search budget, the fits are reproducible up to PyTorch's
## documented non-determinism on multi-threaded CPU operations.
## ============================================================

suppressPackageStartupMessages({
  library(reticulate)
})

## ------------------------------------------------------------
## Lazy module loader. Importing the Python modules at file
## sourcing time would force every user of the library to have
## pycox installed even if they never call the deep fitters.
## We import on first use instead.
## ------------------------------------------------------------
.pycox_state <- new.env(parent = emptyenv())

.ensure_pycox <- function() {
  if (!is.null(.pycox_state$loaded) && .pycox_state$loaded) return(invisible())
  
  .pycox_state$torch       <- import("torch")
  .pycox_state$tt          <- import("torchtuples")
  .pycox_state$pycox       <- import("pycox")
  .pycox_state$pycox_models <- import("pycox.models")
  ## Numpy is imported with convert = FALSE so that the input
  ## chain np$ascontiguousarray(X)$astype("float32") returns a
  ## Python reference rather than auto-converting back to R.
  ## On the OUTPUT path, pycox returns pandas DataFrames that
  ## reticulate auto-converts to R data frames (because pycox
  ## itself was imported with default convert = TRUE), so the
  ## predict closures handle them as plain R objects with no
  ## further conversion needed.
  .pycox_state$np          <- import("numpy", convert = FALSE)
  .pycox_state$loaded      <- TRUE
  invisible()
}


## ------------------------------------------------------------
## Internal helper: convert R training data into the float32
## numpy arrays that pycox expects. pycox is strict about
## dtypes — passing float64 silently degrades training and
## passing the wrong shape raises cryptic Python errors.
## ------------------------------------------------------------
.to_pycox_arrays <- function(X, time, event) {
  np <- .pycox_state$np
  ## With np imported as convert = FALSE (see .ensure_pycox),
  ## np$ascontiguousarray() returns a Python reference rather
  ## than auto-converting back to R, so $astype() chains work.
  ## ascontiguousarray + astype together guarantee a fresh,
  ## contiguous, writable, float32 buffer that survives the
  ## hand-off to torch.from_numpy().
  X_np      <- np$ascontiguousarray(X)$astype("float32")
  durations <- np$ascontiguousarray(as.numeric(time))$astype("float32")
  events    <- np$ascontiguousarray(as.numeric(event))$astype("float32")
  list(X = X_np, durations = durations, events = events)
}


## ------------------------------------------------------------
## Internal helper: split training rows into sub-train and
## validation for the early-stopping signal. We use a stratified
## 80/20 split on the event indicator to keep event proportions
## stable, matching the spirit of the stratified outer CV in
## 04_cv_scaffolding.R.
## ------------------------------------------------------------
.train_val_split <- function(n, event, val_frac = 0.2, seed = 1L) {
  set.seed(seed)
  is_event <- which(event == 1)
  is_cens  <- which(event == 0)
  n_val_event <- max(1, round(length(is_event) * val_frac))
  n_val_cens  <- max(1, round(length(is_cens)  * val_frac))
  val_idx <- c(sample(is_event, n_val_event),
               sample(is_cens,  n_val_cens))
  train_idx <- setdiff(seq_len(n), val_idx)
  list(train = train_idx, val = val_idx)
}


## ------------------------------------------------------------
## Internal helper: build a torchtuples MLPVanilla network with
## the given architecture. Returns a Python network object ready
## to be passed to pycox.models.CoxPH() or CoxTime().
##
## CoxTime expects an extra time input dimension on the first
## layer; the MLPVanillaCoxTime helper handles this. We use the
## right one based on `with_time`.
## ------------------------------------------------------------
.make_network <- function(in_features, num_nodes, dropout, with_time) {
  tt <- .pycox_state$tt
  ## Force num_nodes to a Python list, not a scalar. reticulate
  ## auto-scalarizes length-1 R vectors, so as.integer(c(64)) ends
  ## up as Python int(64). pycox then tries num_nodes[0] and raises
  ## "tuple index out of range" because ints are not subscriptable.
  ## Wrapping in r_to_py(as.list(...)) preserves list semantics for
  ## any length, including 1.
  num_nodes_py <- r_to_py(as.list(as.integer(num_nodes)))
  
  if (with_time) {
    pycox_net <- import("pycox.models.cox_time")
    net <- pycox_net$MLPVanillaCoxTime(
      in_features = as.integer(in_features),
      num_nodes   = num_nodes_py,
      batch_norm  = TRUE,
      dropout     = dropout
    )
  } else {
    ## DeepSurv: same num_nodes handling.
    net <- tt$practical$MLPVanilla(
      in_features  = as.integer(in_features),
      num_nodes    = num_nodes_py,
      out_features = 1L,
      batch_norm   = TRUE,
      dropout      = dropout,
      output_bias  = FALSE
    )
  }
  net
}


## ------------------------------------------------------------
## Internal helper: train one network configuration and return
## (model, val_loss). Used by the random search loop. If
## training raises a Python exception (NaN loss, etc.), returns
## NULL for the model and Inf for the loss so the search can
## skip that configuration.
## ------------------------------------------------------------
.train_one_config <- function(model_class,
                              X_tr, dur_tr, evt_tr,
                              X_val, dur_val, evt_val,
                              num_nodes, dropout, learning_rate,
                              batch_size, n_epochs,
                              with_time) {
  tt <- .pycox_state$tt
  pycox_models <- .pycox_state$pycox_models
  
  net <- .make_network(in_features = ncol(X_tr),
                       num_nodes   = num_nodes,
                       dropout     = dropout,
                       with_time   = with_time)
  
  optimizer <- tt$optim$Adam(lr = learning_rate)
  
  if (model_class == "deepsurv") {
    model <- pycox_models$CoxPH(net, optimizer)
    y_tr_use  <- tuple(dur_tr,  evt_tr)
    y_val_use <- tuple(dur_val, evt_val)
  } else if (model_class == "cox_time") {
    ## CoxTime requires standardized durations. The labtrans
    ## utility handles this and is mandatory for CoxTime.
    ## Pass the durations and events as separate R objects
    ## (the float32 numpy arrays from .to_pycox_arrays), not as
    ## elements of a Python tuple, because indexing a Python
    ## tuple with R's [[i]] uses 0-based semantics inconsistently
    ## across reticulate versions.
    ## NOTE on dtype: labtrans$fit_transform() returns float64
    ## arrays even when fed float32 input. The network expects
    ## float32, so we recast the labtrans output back to float32
    ## via a tiny Python helper. Doing the recast in Python
    ## avoids having to index Python tuples from R, which has
    ## been the source of repeated bugs in this wrapper.
    labtrans <- pycox_models$CoxTime$label_transform()
    y_tr_raw  <- labtrans$fit_transform(dur_tr,  evt_tr)
    y_val_raw <- labtrans$transform(dur_val, evt_val)
    py_recast <- reticulate::py_run_string(
      "def _recast_yt(y):
    import numpy as np
    return (np.ascontiguousarray(y[0]).astype('float32'),
            np.ascontiguousarray(y[1]).astype('float32'))",
      convert = FALSE
    )
    y_tr_t   <- py_recast$`_recast_yt`(y_tr_raw)
    y_val_t  <- py_recast$`_recast_yt`(y_val_raw)
    model    <- pycox_models$CoxTime(net, optimizer, labtrans = labtrans)
    y_tr_use  <- y_tr_t
    y_val_use <- y_val_t
  } else {
    stop("internal error: unknown model_class")
  }
  
  callbacks <- list(tt$callbacks$EarlyStopping(patience = 10L))
  
  log <- tryCatch(
    model$fit(
      X_tr, y_tr_use,
      batch_size = as.integer(batch_size),
      epochs     = as.integer(n_epochs),
      callbacks  = callbacks,
      val_data   = tuple(X_val, y_val_use),
      verbose    = FALSE
    ),
    error = function(e) {
      message(sprintf("  pycox fit failed for one config: %s",
                      conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(log)) {
    return(list(model = NULL, val_loss = Inf))
  }
  
  ## Extract final validation loss from the training log.
  val_loss <- tryCatch(
    {
      to_pandas <- log$to_pandas()
      tail(to_pandas$val_loss, 1)
    },
    error = function(e) Inf
  )
  if (!is.finite(val_loss)) val_loss <- Inf
  
  list(model = model, val_loss = val_loss)
}


## ============================================================
## fit_deepsurv()
##
## DeepSurv via pycox.models.CoxPH. Section III-B.
##
## Random hyperparameter search over:
##   - learning rate     {1e-4, 1e-3, 1e-2}
##   - dropout           {0.0, 0.1, 0.3, 0.5}
##   - hidden width      {32, 64, 128}
##   - number of layers  {1, 2, 3}
##
## Default search budget n_random_search = 8 configurations.
## Each is trained with early stopping on a 20% validation split.
## The best-validation-loss configuration is retrained on the
## full training fold and used for prediction.
## ============================================================
fit_deepsurv <- function(X_train, time_train, event_train,
                         inner_cv = NULL, seed = NULL,
                         n_random_search = 8,
                         n_epochs = 200,
                         batch_size = 64L, ...) {
  t0 <- Sys.time()
  .ensure_pycox()
  if (!is.null(seed)) {
    set.seed(seed)
    .pycox_state$torch$manual_seed(as.integer(seed))
  }
  
  ## --- Train/validation split for early stopping ---
  split <- .train_val_split(nrow(X_train), event_train,
                            val_frac = 0.2, seed = seed %||% 1L)
  arr_tr  <- .to_pycox_arrays(X_train[split$train, , drop = FALSE],
                              time_train[split$train],
                              event_train[split$train])
  arr_val <- .to_pycox_arrays(X_train[split$val, , drop = FALSE],
                              time_train[split$val],
                              event_train[split$val])
  X_tr    <- arr_tr$X
  X_val   <- arr_val$X
  dur_tr  <- arr_tr$durations
  evt_tr  <- arr_tr$events
  dur_val <- arr_val$durations
  evt_val <- arr_val$events
  
  ## --- Build random search grid ---
  lr_grid       <- c(1e-4, 1e-3, 1e-2)
  dropout_grid  <- c(0.0, 0.1, 0.3, 0.5)
  width_grid    <- c(32L, 64L, 128L)
  depth_grid    <- c(1L, 2L, 3L)
  
  configs <- vector("list", n_random_search)
  for (s in seq_len(n_random_search)) {
    width <- sample(width_grid, 1)
    depth <- sample(depth_grid, 1)
    configs[[s]] <- list(
      learning_rate = sample(lr_grid, 1),
      dropout       = sample(dropout_grid, 1),
      num_nodes     = rep(width, depth)
    )
  }
  
  ## --- Train each config ---
  best_loss <- Inf
  best_config <- NULL
  for (s in seq_along(configs)) {
    cfg <- configs[[s]]
    res <- .train_one_config(
      model_class   = "deepsurv",
      X_tr   = X_tr,   dur_tr  = dur_tr,  evt_tr  = evt_tr,
      X_val  = X_val,  dur_val = dur_val, evt_val = evt_val,
      num_nodes     = cfg$num_nodes,
      dropout       = cfg$dropout,
      learning_rate = cfg$learning_rate,
      batch_size    = batch_size,
      n_epochs      = n_epochs,
      with_time     = FALSE
    )
    if (res$val_loss < best_loss) {
      best_loss <- res$val_loss
      best_config <- cfg
    }
  }
  
  if (is.null(best_config)) {
    stop("fit_deepsurv: all random search configurations failed")
  }
  
  ## --- Refit on the full training fold with the best config ---
  arr_full <- .to_pycox_arrays(X_train, time_train, event_train)
  X_full <- arr_full$X
  ## y_full is built as a fresh Python tuple right where it is
  ## consumed; we no longer carry tuples around as R-side state
  ## because indexing them back open is fragile across reticulate
  ## versions.
  y_full <- tuple(arr_full$durations, arr_full$events)
  
  net <- .make_network(in_features = ncol(X_full),
                       num_nodes   = best_config$num_nodes,
                       dropout     = best_config$dropout,
                       with_time   = FALSE)
  optimizer <- .pycox_state$tt$optim$Adam(lr = best_config$learning_rate)
  raw_fit <- .pycox_state$pycox_models$CoxPH(net, optimizer)
  raw_fit$fit(
    X_full, y_full,
    batch_size = as.integer(batch_size),
    epochs     = as.integer(n_epochs),
    callbacks  = list(.pycox_state$tt$callbacks$EarlyStopping(patience = 10L)),
    val_data   = tuple(X_val, tuple(dur_val, evt_val)),
    verbose    = FALSE
  )
  ## CoxPH needs the baseline hazards computed on the training
  ## data after fitting; pycox provides a convenience method.
  raw_fit$compute_baseline_hazards()
  
  np <- .pycox_state$np
  
  predict_risk <- function(X_new) {
    X_new_np <- .pycox_state$np$ascontiguousarray(X_new)$astype("float32")
    ## predict() returns an n_new x 1 R matrix because the
    ## auto-converter handles the numpy return cleanly.
    pred <- raw_fit$predict(X_new_np)
    as.numeric(pred)
  }
  
  predict_survival <- function(X_new, times) {
    X_new_np <- .pycox_state$np$ascontiguousarray(X_new)$astype("float32")
    ## predict_surv_df returns a pandas DataFrame indexed by
    ## time, with one column per subject. We evaluate it at
    ## the requested times.
    ## predict_surv_df() returns an R data frame thanks to
    ## reticulate auto-conversion on the pycox import. Row names
    ## are the time grid (in days), column names are subject
    ## indices, and cells are predicted survival probabilities.
    surv_df <- raw_fit$predict_surv_df(X_new_np)
    surv_mat <- as.matrix(surv_df)
    base_times <- as.numeric(rownames(surv_df))
    
    ## Step-function interpolation onto the requested grid.
    ## Two edge cases that matter for the benchmark:
    ##   - times[k] BEFORE pycox's grid starts: use the earliest
    ##     available row of surv_mat rather than constant 1.
    ##     Returning constant 1 destroys all between-subject
    ##     variation and breaks the calibration regression.
    ##   - times[k] AFTER pycox's grid ends: use the last
    ##     available row rather than 0. Returning 0 inflates
    ##     IBS contributions because the squared error becomes
    ##     (1 - 0)^2 = 1 for every still-at-risk subject.
    n_new <- nrow(X_new)
    n_grid <- length(base_times)
    out <- matrix(1, nrow = n_new, ncol = length(times))
    for (k in seq_along(times)) {
      idx <- findInterval(times[k], base_times)
      idx <- max(1L, min(idx, n_grid))   # clamp to [1, n_grid]
      out[, k] <- surv_mat[idx, ]
    }
    out
  }
  
  list(
    model_class      = "deepsurv",
    raw_fit          = raw_fit,
    coef             = NULL,
    predict_risk     = predict_risk,
    predict_survival = predict_survival,
    hyperparams      = c(best_config,
                         list(n_random_search = n_random_search,
                              n_epochs        = n_epochs,
                              batch_size      = batch_size,
                              best_val_loss   = best_loss)),
    fit_time_sec     = as.numeric(Sys.time() - t0, units = "secs")
  )
}


## ============================================================
## fit_cox_time()
##
## Cox-Time via pycox.models.CoxTime. Section III-B. Same
## random search structure as DeepSurv. The internal handling
## of the time input is delegated to pycox's MLPVanillaCoxTime
## network, and the durations are standardized via the labtrans
## utility (mandatory for CoxTime per the pycox documentation).
## ============================================================
fit_cox_time <- function(X_train, time_train, event_train,
                         inner_cv = NULL, seed = NULL,
                         n_random_search = 8,
                         n_epochs = 200,
                         batch_size = 64L, ...) {
  t0 <- Sys.time()
  .ensure_pycox()
  if (!is.null(seed)) {
    set.seed(seed)
    .pycox_state$torch$manual_seed(as.integer(seed))
  }
  
  split <- .train_val_split(nrow(X_train), event_train,
                            val_frac = 0.2, seed = seed %||% 1L)
  arr_tr  <- .to_pycox_arrays(X_train[split$train, , drop = FALSE],
                              time_train[split$train],
                              event_train[split$train])
  arr_val <- .to_pycox_arrays(X_train[split$val, , drop = FALSE],
                              time_train[split$val],
                              event_train[split$val])
  X_tr    <- arr_tr$X
  X_val   <- arr_val$X
  dur_tr  <- arr_tr$durations
  evt_tr  <- arr_tr$events
  dur_val <- arr_val$durations
  evt_val <- arr_val$events
  
  lr_grid       <- c(1e-4, 1e-3, 1e-2)
  dropout_grid  <- c(0.0, 0.1, 0.3, 0.5)
  width_grid    <- c(32L, 64L, 128L)
  depth_grid    <- c(1L, 2L, 3L)
  
  configs <- vector("list", n_random_search)
  for (s in seq_len(n_random_search)) {
    width <- sample(width_grid, 1)
    depth <- sample(depth_grid, 1)
    configs[[s]] <- list(
      learning_rate = sample(lr_grid, 1),
      dropout       = sample(dropout_grid, 1),
      num_nodes     = rep(width, depth)
    )
  }
  
  best_loss <- Inf
  best_config <- NULL
  for (s in seq_along(configs)) {
    cfg <- configs[[s]]
    res <- .train_one_config(
      model_class   = "cox_time",
      X_tr   = X_tr,   dur_tr  = dur_tr,  evt_tr  = evt_tr,
      X_val  = X_val,  dur_val = dur_val, evt_val = evt_val,
      num_nodes     = cfg$num_nodes,
      dropout       = cfg$dropout,
      learning_rate = cfg$learning_rate,
      batch_size    = batch_size,
      n_epochs      = n_epochs,
      with_time     = TRUE
    )
    if (res$val_loss < best_loss) {
      best_loss <- res$val_loss
      best_config <- cfg
    }
  }
  
  if (is.null(best_config)) {
    stop("fit_cox_time: all random search configurations failed")
  }
  
  ## Refit on full training fold with best config.
  arr_full <- .to_pycox_arrays(X_train, time_train, event_train)
  X_full <- arr_full$X
  
  net <- .make_network(in_features = ncol(X_full),
                       num_nodes   = best_config$num_nodes,
                       dropout     = best_config$dropout,
                       with_time   = TRUE)
  optimizer <- .pycox_state$tt$optim$Adam(lr = best_config$learning_rate)
  
  labtrans <- .pycox_state$pycox_models$CoxTime$label_transform()
  y_full_raw <- labtrans$fit_transform(arr_full$durations, arr_full$events)
  y_val_raw  <- labtrans$transform(arr_val$durations, arr_val$events)
  ## Recast labtrans output back to float32 (see notes in
  ## .train_one_config). Same Python helper.
  py_recast <- reticulate::py_run_string(
    "def _recast_yt(y):
    import numpy as np
    return (np.ascontiguousarray(y[0]).astype('float32'),
            np.ascontiguousarray(y[1]).astype('float32'))",
    convert = FALSE
  )
  y_full_t <- py_recast$`_recast_yt`(y_full_raw)
  y_val_t  <- py_recast$`_recast_yt`(y_val_raw)
  
  raw_fit <- .pycox_state$pycox_models$CoxTime(net, optimizer,
                                               labtrans = labtrans)
  raw_fit$fit(
    X_full, y_full_t,
    batch_size = as.integer(batch_size),
    epochs     = as.integer(n_epochs),
    callbacks  = list(.pycox_state$tt$callbacks$EarlyStopping(patience = 10L)),
    val_data   = tuple(X_val, y_val_t),
    verbose    = FALSE
  )
  raw_fit$compute_baseline_hazards()
  
  np <- .pycox_state$np
  
  predict_risk <- function(X_new) {
    X_new_np <- .pycox_state$np$ascontiguousarray(X_new)$astype("float32")
    ## For Cox-Time, "risk" is a model-internal quantity that
    ## depends on time. We summarize it as the time-integrated
    ## log relative risk at the median training event time, so
    ## the result is a single number per subject usable by the
    ## C-index. This is the standard approach for Cox-Time
    ## risk ranking (Kvamme et al. 2019, Section 4.1).
    surv_df <- raw_fit$predict_surv_df(X_new_np)
    surv_mat <- as.matrix(surv_df)
    ## Use -log(S(t_med)) as the risk score; higher means higher
    ## hazard, matching the harrell_c sign convention.
    n_t <- nrow(surv_mat)
    t_med_row <- max(1, floor(n_t / 2))
    surv_med <- surv_mat[t_med_row, ]
    pmax(-log(pmax(surv_med, 1e-8)), 0)
  }
  
  predict_survival <- function(X_new, times) {
    X_new_np <- .pycox_state$np$ascontiguousarray(X_new)$astype("float32")
    surv_df <- raw_fit$predict_surv_df(X_new_np)
    surv_mat <- as.matrix(surv_df)
    base_times <- as.numeric(rownames(surv_df))
    
    n_new <- nrow(X_new)
    n_grid <- length(base_times)
    out <- matrix(1, nrow = n_new, ncol = length(times))
    for (k in seq_along(times)) {
      idx <- findInterval(times[k], base_times)
      idx <- max(1L, min(idx, n_grid))   # clamp to [1, n_grid]
      out[, k] <- surv_mat[idx, ]
    }
    out
  }
  
  list(
    model_class      = "cox_time",
    raw_fit          = raw_fit,
    coef             = NULL,
    predict_risk     = predict_risk,
    predict_survival = predict_survival,
    hyperparams      = c(best_config,
                         list(n_random_search = n_random_search,
                              n_epochs        = n_epochs,
                              batch_size      = batch_size,
                              best_val_loss   = best_loss)),
    fit_time_sec     = as.numeric(Sys.time() - t0, units = "secs")
  )
}


## ------------------------------------------------------------
## %||% operator for "use right side if left is NULL". Defined
## here for self-containment in case the user has not loaded
## rlang or purrr.
## ------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a