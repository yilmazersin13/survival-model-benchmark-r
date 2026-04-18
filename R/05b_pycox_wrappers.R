## ============================================================
## 05b_pycox_wrappers.R



suppressPackageStartupMessages({
  library(reticulate)
})


.pycox_state <- new.env(parent = emptyenv())

.ensure_pycox <- function() {
  if (!is.null(.pycox_state$loaded) && .pycox_state$loaded) return(invisible())
  
  .pycox_state$torch       <- import("torch")
  .pycox_state$tt          <- import("torchtuples")
  .pycox_state$pycox       <- import("pycox")
  .pycox_state$pycox_models <- import("pycox.models")

  .pycox_state$np          <- import("numpy", convert = FALSE)
  .pycox_state$loaded      <- TRUE
  invisible()
}


.to_pycox_arrays <- function(X, time, event) {
  np <- .pycox_state$np

  X_np      <- np$ascontiguousarray(X)$astype("float32")
  durations <- np$ascontiguousarray(as.numeric(time))$astype("float32")
  events    <- np$ascontiguousarray(as.numeric(event))$astype("float32")
  list(X = X_np, durations = durations, events = events)
}


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



.make_network <- function(in_features, num_nodes, dropout, with_time) {
  tt <- .pycox_state$tt

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

    surv_df <- raw_fit$predict_surv_df(X_new_np)
    surv_mat <- as.matrix(surv_df)

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



`%||%` <- function(a, b) if (is.null(a)) b else a
