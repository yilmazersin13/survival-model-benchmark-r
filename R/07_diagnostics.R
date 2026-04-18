## ============================================================
## 07_diagnostics.R
##


suppressPackageStartupMessages({
  library(survival)
})

##
## Arguments:
##   fit       a fitted-model object from 05_model_fitters.R
##             (or 05b_ben_cox.R)
##   p_raw     the post-preprocessing predictor count, which is
##             the number of columns in the X_train matrix that
##             was passed to the fitter. Used as the fallback
##             value for models without a sparsity notion.
##
## Returns a numeric scalar p_eff. Always positive; floored at

compute_p_eff <- function(fit, p_raw) {
  if (is.null(fit$model_class)) {
    stop("compute_p_eff: fit object has no $model_class slot")
  }
  if (!is.numeric(p_raw) || p_raw < 1) {
    stop("compute_p_eff: p_raw must be a positive number")
  }
  
  cls <- fit$model_class
  
  p_eff <- switch(
    cls,
    "cox"      = p_raw,
    "rsf"      = p_raw,
    "deepsurv" = p_raw,
    "cox_time" = p_raw,
    
    "lasso_cox" = {
      ## Number of nonzero coefficients at the chosen lambda.
      ## The fitter stores the coefficient vector at lambda.min
      ## in $coef.
      if (is.null(fit$coef)) {
        stop("compute_p_eff: lasso_cox fit missing $coef")
      }
      sum(fit$coef != 0)
    },
    
    "ben_cox" = {

      if (is.null(fit$beta) || is.null(fit$beta$quantiles)) {
        stop("compute_p_eff: ben_cox fit missing $beta$quantiles")
      }
      q <- fit$beta$quantiles
      sum(q[1, ] > 0 | q[5, ] < 0)
    },
    
    stop(sprintf("compute_p_eff: unknown model_class '%s'", cls))
  )
  
  ## Floor at 1 so that EPV is always finite. A model that
  ## shrinks every coefficient to exactly zero is degenerate
  ## but should not crash the diagnostics pipeline.
  max(p_eff, 1)
}



compute_epv <- function(event_train_aug, fit, p_raw) {
  if (!all(event_train_aug %in% c(0L, 1L, 0, 1))) {
    stop("compute_epv: event_train_aug must be 0/1")
  }
  n_events <- sum(event_train_aug == 1)
  p_eff <- compute_p_eff(fit, p_raw)
  
  list(
    epv      = n_events / p_eff,
    epv_raw  = n_events / p_raw,
    n_events = as.integer(n_events),
    p_eff    = as.integer(p_eff),
    p_raw    = as.integer(p_raw)
  )
}

## Arguments:
##   X_train      numeric matrix, n_train x p_raw, the
##                preprocessed training-fold design matrix
##   time_train   numeric vector of training-fold times
##   event_train  integer vector of training-fold events
##
## Returns a list with:
##   $global_p    numeric scalar in (0,1), or NA if the test
##                could not be computed
##   $chisq       numeric scalar, the global chi-squared statistic
##   $df          integer, degrees of freedom of the global test
##   $converged   logical, TRUE if the unpenalized Cox fit
##                converged without warnings
##   $events_per_p  numeric, n_events / p_raw, recorded for
##                  downstream flagging of low-power cells

ph_test_global <- function(X_train, time_train, event_train) {
  if (nrow(X_train) != length(time_train)) {
    stop("ph_test_global: X_train and time_train length mismatch")
  }
  if (length(time_train) != length(event_train)) {
    stop("ph_test_global: time_train and event_train length mismatch")
  }
  
  p_raw <- ncol(X_train)
  n_events <- sum(event_train == 1)
  events_per_p <- n_events / p_raw
  
  na_result <- list(
    global_p     = NA_real_,
    chisq        = NA_real_,
    df           = NA_integer_,
    converged    = FALSE,
    events_per_p = events_per_p
  )
  

  if (is.null(colnames(X_train))) {
    colnames(X_train) <- paste0("V", seq_len(p_raw))
  }
  df <- data.frame(time = time_train, event = event_train,
                   as.data.frame(X_train),
                   check.names = FALSE)
  
  formula_str <- sprintf("Surv(time, event) ~ %s",
                         paste0("`", colnames(X_train), "`",
                                collapse = " + "))
  
  ## Fit unpenalized Cox, catching all warnings and errors.
  ## A non-converged fit is treated as a failure and we return
  ## the NA result.
  fit_status <- list(converged = TRUE)
  cox_fit <- withCallingHandlers(
    tryCatch(
      coxph(as.formula(formula_str), data = df, model = FALSE,
            x = FALSE, y = FALSE),
      error = function(e) NULL
    ),
    warning = function(w) {
      fit_status$converged <<- FALSE
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(cox_fit) || !fit_status$converged) {
    return(na_result)
  }
  if (!isTRUE(cox_fit$iter > 0)) {
    return(na_result)
  }
  
  ## Run the Grambsch-Therneau test. cox.zph returns a list with
  ## a `table` matrix whose last row is the global test.
  zph <- tryCatch(
    cox.zph(cox_fit, transform = "km", global = TRUE),
    error = function(e) NULL
  )
  if (is.null(zph)) return(na_result)
  
  global_row <- zph$table[nrow(zph$table), , drop = TRUE]
  if (any(!is.finite(global_row))) return(na_result)
  
  list(
    global_p     = unname(global_row["p"]),
    chisq        = unname(global_row["chisq"]),
    df           = as.integer(unname(global_row["df"])),
    converged    = TRUE,
    events_per_p = events_per_p
  )
}
