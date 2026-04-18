
## 02_preprocessing.R
##
## Fold-aware preprocessing implementing Section IV-B of the
## paper. The contract is strict:
##

## Public functions:
##   fit_preprocessor()    learn statistics from a training subset
##   apply_preprocessor()  transform a row subset using stored stats
##   preprocess_fold()     convenience wrapper: fit on train_idx,
##                         apply to both train_idx and test_idx,
##                         return the two transformed matrices

.statistical_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  tab <- table(x)
  names(tab)[which.max(tab)]
}


fit_preprocessor <- function(dataset, train_idx,
                             apply_tcga_filter = TRUE) {
  if (!all(c("X", "feature_names", "feature_types") %in% names(dataset))) {
    stop("fit_preprocessor: input does not look like a loader object")
  }
  
  X_train <- dataset$X[train_idx, , drop = FALSE]
  feature_names <- dataset$feature_names
  feature_types <- dataset$feature_types
  cat_levels_attr <- attr(dataset$X, "categorical_levels")
  
  cont_idx_all <- which(feature_types == "continuous")
  cat_idx      <- which(feature_types == "categorical")
  

  if (apply_tcga_filter &&
      identical(dataset$dataset_name, "TCGA-BRCA") &&
      length(cont_idx_all) > 500) {
    
    cont_block <- X_train[, cont_idx_all, drop = FALSE]
    train_vars <- apply(cont_block, 2, function(col) {
      var(col, na.rm = TRUE)
    })
    train_vars[is.na(train_vars)] <- -Inf
    keep_local <- order(train_vars, decreasing = TRUE)[1:500]
    cont_idx <- cont_idx_all[sort(keep_local)]
    message(sprintf(
      "fit_preprocessor: TCGA-BRCA filter applied, retained %d of %d continuous probes",
      length(cont_idx), length(cont_idx_all)))
  } else {
    cont_idx <- cont_idx_all
  }
  
  ## --- Continuous statistics (training fold only) ---
  if (length(cont_idx) > 0) {
    cont_block <- X_train[, cont_idx, drop = FALSE]
    cont_means   <- apply(cont_block, 2, mean,   na.rm = TRUE)
    cont_sds     <- apply(cont_block, 2, sd,     na.rm = TRUE)
    cont_medians <- apply(cont_block, 2, median, na.rm = TRUE)

    cont_sds[!is.finite(cont_sds) | cont_sds == 0] <- 1
    

    cont_medians[!is.finite(cont_medians)] <- 0
  } else {
    cont_means   <- numeric(0)
    cont_sds     <- numeric(0)
    cont_medians <- numeric(0)
  }
  
  ## --- Categorical statistics (training fold only) ---
  cat_levels <- list()
  cat_modes  <- character(0)
  if (length(cat_idx) > 0) {
    for (j in cat_idx) {
      col_name <- feature_names[j]

      lvl_vec <- cat_levels_attr[[col_name]]
      if (is.null(lvl_vec)) {

        train_codes <- X_train[, j]
        train_chars <- as.character(train_codes)
      } else {
        train_codes <- X_train[, j]
        train_chars <- ifelse(is.na(train_codes), NA_character_,
                              lvl_vec[train_codes])
      }
      
      observed <- unique(train_chars[!is.na(train_chars)])
      cat_levels[[col_name]] <- observed
      cat_modes[col_name]    <- .statistical_mode(train_chars)
    }
  }
  

  cont_names <- feature_names[cont_idx]
  onehot_names <- character(0)
  for (j in cat_idx) {
    col_name <- feature_names[j]
    lvls <- cat_levels[[col_name]]
    if (length(lvls) > 0) {
      onehot_names <- c(onehot_names, paste0(col_name, "__", lvls))
    }
  }
  output_colnames <- c(cont_names, onehot_names)
  
  list(
    continuous_idx  = cont_idx,
    categorical_idx = cat_idx,
    cont_means      = cont_means,
    cont_sds        = cont_sds,
    cont_medians    = cont_medians,
    cat_levels      = cat_levels,
    cat_modes       = cat_modes,
    output_colnames = output_colnames,
    dataset_name    = dataset$dataset_name,
    feature_names   = feature_names,
    cat_levels_attr = cat_levels_attr
  )
}



apply_preprocessor <- function(prep, dataset, row_idx) {
  X_sub <- dataset$X[row_idx, , drop = FALSE]
  n_out <- length(row_idx)
  
  ## --- Continuous block ---
  if (length(prep$continuous_idx) > 0) {
    cont_block <- X_sub[, prep$continuous_idx, drop = FALSE]
    ## Impute with training-fold medians.
    for (k in seq_len(ncol(cont_block))) {
      col <- cont_block[, k]
      na_mask <- is.na(col)
      if (any(na_mask)) {
        col[na_mask] <- prep$cont_medians[k]
        cont_block[, k] <- col
      }
    }
    ## Standardize with training-fold means and SDs.
    cont_block <- sweep(cont_block, 2, prep$cont_means, "-")
    cont_block <- sweep(cont_block, 2, prep$cont_sds,   "/")
  } else {
    cont_block <- matrix(numeric(0), nrow = n_out, ncol = 0)
  }

  n_onehot <- length(prep$output_colnames) - length(prep$continuous_idx)
  if (n_onehot > 0) {
    onehot_block <- matrix(0, nrow = n_out, ncol = n_onehot)
    out_col <- 1L
    for (j in prep$categorical_idx) {
      col_name <- prep$feature_names[j]
      lvls     <- prep$cat_levels[[col_name]]
      if (length(lvls) == 0) next
      
      lvl_vec <- prep$cat_levels_attr[[col_name]]
      codes <- X_sub[, j]
      if (is.null(lvl_vec)) {
        chars <- as.character(codes)
      } else {
        chars <- ifelse(is.na(codes), NA_character_, lvl_vec[codes])
      }
      ## Impute missing categorical values with the training mode.
      chars[is.na(chars)] <- prep$cat_modes[col_name]
      
      for (lv in lvls) {
        onehot_block[, out_col] <- as.integer(chars == lv)
        out_col <- out_col + 1L
      }
    }
  } else {
    onehot_block <- matrix(numeric(0), nrow = n_out, ncol = 0)
  }
  
  out <- cbind(cont_block, onehot_block)
  colnames(out) <- prep$output_colnames
  out
}


preprocess_fold <- function(dataset, train_idx, test_idx, ...) {
  prep <- fit_preprocessor(dataset, train_idx, ...)
  X_train <- apply_preprocessor(prep, dataset, train_idx)
  X_test  <- apply_preprocessor(prep, dataset, test_idx)
  
  list(
    X_train     = X_train,
    X_test      = X_test,
    time_train  = dataset$time[train_idx],
    event_train = as.integer(dataset$event[train_idx]),
    time_test   = dataset$time[test_idx],
    event_test  = as.integer(dataset$event[test_idx]),
    prep        = prep
  )
}
