## ============================================================
## 02_preprocessing.R
##
## Fold-aware preprocessing implementing Section IV-B of the
## paper. The contract is strict:
##
##   - Continuous covariates are standardized to zero mean and
##     unit variance using statistics computed on TRAINING rows
##     only.
##   - Categorical covariates are one-hot encoded, with the set
##     of levels determined from TRAINING rows only. Levels that
##     appear in the test fold but not in training are mapped to
##     all-zero rows (the "unseen level" convention).
##   - Missing continuous values are imputed with the
##     TRAINING-fold median.
##   - Missing categorical values are imputed with the
##     TRAINING-fold mode (the most frequent observed level).
##   - The TCGA-BRCA top-500-variance probe filter from
##     Section IV-B, when applicable, is computed on training
##     rows only and applied to test rows.
##
## The two-stage interface mirrors how scikit-learn or the R
## `recipes` package handles this: a `fit` step that learns
## statistics from training data and stores them, and an `apply`
## step that uses the stored statistics to transform any new
## data (training fold itself, test fold, or future data).
##
## This separation is what guarantees no leakage: the statistics
## are computed once, on training rows only, and then applied to
## both training and test. There is no path by which test
## information can enter the statistics.
##
## Public functions:
##   fit_preprocessor()    learn statistics from a training subset
##   apply_preprocessor()  transform a row subset using stored stats
##   preprocess_fold()     convenience wrapper: fit on train_idx,
##                         apply to both train_idx and test_idx,
##                         return the two transformed matrices
## ============================================================

## ------------------------------------------------------------
## Internal helper: compute the statistical mode of a categorical
## vector, ignoring NA. Ties are broken by the order in which
## levels first appear, which is deterministic given a fixed
## input ordering.
## ------------------------------------------------------------
.statistical_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  tab <- table(x)
  names(tab)[which.max(tab)]
}


## ============================================================
## fit_preprocessor()
##
## Learn preprocessing statistics from a training subset of a
## dataset object. Returns a "preprocessor" object that records
## every quantity needed to transform new data, plus enough
## metadata for downstream code to know what columns to expect.
##
## Arguments:
##   dataset       a dataset object from 01_data_loaders.R
##                 (or its augmented counterpart from
##                 03_censoring_augmentation.R)
##   train_idx     integer vector, the row indices of the
##                 training fold (positions in the original
##                 n-row dataset, as produced by
##                 04_cv_scaffolding.R)
##   apply_tcga_filter  logical; if TRUE and the dataset is
##                 TCGA-BRCA with at least 500 continuous
##                 covariates, retain only the 500 most variable
##                 (computed on the training fold). Default
##                 follows the spirit of Section IV-B and is
##                 controlled by the caller.
##
## Returns a list with:
##   $continuous_idx    integer vector of original column indices
##                      that are continuous (and survived the
##                      TCGA filter, if applied)
##   $categorical_idx   integer vector of original column indices
##                      that are categorical
##   $cont_means        numeric vector of training-fold means,
##                      one per continuous covariate (post-filter)
##   $cont_sds          numeric vector of training-fold SDs,
##                      one per continuous covariate (post-filter)
##   $cont_medians      numeric vector of training-fold medians
##                      for imputation (post-filter)
##   $cat_levels        list keyed by categorical column name,
##                      each element a character vector of levels
##                      observed in the training fold
##   $cat_modes         named character vector of training-fold
##                      modes for imputation
##   $output_colnames   character vector of the columns that
##                      apply_preprocessor() will produce, in
##                      order. Continuous first, then one-hot
##                      indicators in the order
##                      <feature>__<level>.
##   $dataset_name      passed through for diagnostic messages
## ============================================================
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
  
  ## --- Optional TCGA-BRCA top-500-variance filter ---
  ## Section IV-B: "the number of molecular features is reduced
  ## to the p = 500 most variable probes prior to modeling."
  ## We compute the variances on the training fold only, never
  ## on the full dataset, so that the filter respects fold
  ## isolation. The filter only fires when the dataset is
  ## TCGA-BRCA and the continuous block has more than 500 columns,
  ## so it is a no-op for METABRIC, SUPPORT, and the clinical-only
  ## TCGA-BRCA cohort the user currently has loaded.
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
    
    ## Guard against zero-variance columns: replace SD = 0 with 1
    ## so that the standardization step does not divide by zero.
    ## A zero-variance column carries no information anyway, but
    ## we keep it rather than dropping it so that the column
    ## structure is identical across folds.
    cont_sds[!is.finite(cont_sds) | cont_sds == 0] <- 1
    
    ## Guard against all-NA continuous columns: replace NA median
    ## with 0, which is the post-standardization mean. The same
    ## logic applies as above.
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
      ## Recover the original character labels from the integer
      ## codes produced by the loader. The loader stored the level
      ## vector as an attribute of the X matrix.
      lvl_vec <- cat_levels_attr[[col_name]]
      if (is.null(lvl_vec)) {
        ## Fallback: column is numeric integer-coded with no
        ## attached levels. Treat the integers themselves as
        ## levels.
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
  
  ## --- Build the output column name vector ---
  ## Continuous columns keep their original feature names. Each
  ## categorical column expands into one indicator per training
  ## level, named "<feature>__<level>". This is the order in which
  ## apply_preprocessor() will write columns into the output matrix.
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


## ============================================================
## apply_preprocessor()
##
## Transform a row subset of a dataset using a preprocessor
## fitted on (typically) a different subset. The output matrix
## has columns in exactly the order recorded in
## prep$output_colnames, regardless of which rows are passed in.
##
## Arguments:
##   prep      preprocessor object from fit_preprocessor()
##   dataset   the same dataset object the preprocessor was fit on
##             (or any object with the same column structure)
##   row_idx   integer vector of rows to transform; positions in
##             the original n-row dataset
##
## Returns a numeric matrix with length(row_idx) rows and
## length(prep$output_colnames) columns. Suitable as direct
## input to glmnet, randomForestSRC, etc.
## ============================================================
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
  
  ## --- One-hot block ---
  ## Build column-by-column using the levels recorded at fit time.
  ## Levels present in row_idx but not in training are mapped to
  ## all-zero rows for that feature (the "unseen level"
  ## convention). Missing values are imputed with the training-fold
  ## mode and then encoded normally.
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


## ============================================================
## preprocess_fold()
##
## Convenience wrapper that ties fit_preprocessor() and
## apply_preprocessor() together for a single outer fold. This
## is the function the main loop in main.R will actually call.
##
## Arguments:
##   dataset    a dataset object (possibly already augmented by
##              augment_dataset() from 03_censoring_augmentation.R)
##   train_idx  integer vector, training rows in the original
##              dataset
##   test_idx   integer vector, test rows in the original dataset
##   ...        passed through to fit_preprocessor()
##
## Returns a list with:
##   $X_train       n_train x p_eff_raw matrix
##   $X_test        n_test x p_eff_raw matrix
##   $time_train    numeric vector of training times
##   $event_train   integer vector of training events
##   $time_test     numeric vector of test times
##   $event_test    integer vector of test events
##   $prep          the preprocessor object (so downstream code
##                  can read $output_colnames, etc.)
##
## All four (time, event) vectors are taken from the dataset as
## passed in, so if the dataset has been augmented, the augmented
## times and events are what you get. This is the intended
## composition: augment first (on the full dataset, with a fold-
## specific seed), then preprocess on the augmented object.
## ============================================================
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