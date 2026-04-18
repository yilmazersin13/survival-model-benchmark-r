## ============================================================
## 07_diagnostics.R
##
## Three diagnostic computations called out in the paper:
##
##   1. Effective covariate dimension  p_eff   (Section III-C)
##   2. Events-per-variable ratio      EPV     (Section IV-E,
##                                              Equation (5))
##   3. Grambsch-Therneau global PH test       (Section III-C)
##
## All three are model-aware in the sense that they need to know
## which model produced a fit. The convention in this library is
## that every fitter from 05_model_fitters.R returns an object
## with a `$model_class` slot containing one of:
##
##   "cox"        unpenalized Cox
##   "lasso_cox"  lasso-penalized Cox
##   "ben_cox"    Bayesian elastic net Cox
##   "rsf"        random survival forest
##   "deepsurv"   DeepSurv
##   "cox_time"   Cox-Time
##
## The diagnostics dispatch on this slot. The fitters in Step 7
## will be written to populate it.
## ============================================================

suppressPackageStartupMessages({
  library(survival)
})


## ============================================================
## compute_p_eff()
##
## Effective covariate dimension as defined in Section III-C.
## The definition is model-specific:
##
##   - Unpenalized Cox, RSF, DeepSurv, Cox-Time:
##       p_eff = p_raw  (the post-preprocessing predictor count)
##     This is a conservative proxy, not a formal degrees-of-
##     freedom estimate. Section III-C acknowledges the asymmetry
##     and Section VII-D will flag it as a limitation.
##
##   - Lasso-Cox:
##       p_eff = number of nonzero coefficients at the cross-
##       validated penalty lambda. This is the standard lasso
##       degrees-of-freedom estimate (Zou, Hastie & Tibshirani,
##       2007).
##
##   - BEN-Cox:
##       p_eff = posterior expected number of "unshrunk"
##       coefficients. We define a coefficient as effectively
##       unshrunk when its 95 percent credible interval excludes
##       zero, and average this over posterior draws. This gives
##       a Bayesian analog of the lasso's nonzero count and is
##       the estimate referred to in Section III-C.
##
##       The alternative shrinkage-factor formulation
##         p_eff = E[ sum_j (1 - kappa_j) | data ]
##       requires kappa_j to be defined explicitly from the prior
##       and posterior variance, which is awkward for the elastic
##       net construction in 05b_ben_cox.R because the prior has
##       two scale parameters (lambda1 and lambda2) entering the
##       coefficient variance jointly. We therefore use the
##       credible-interval rule, which is interpretable, robust,
##       and what the user's existing BEN-Cox code already
##       computes when reporting "significant genes."
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
## 1 to avoid divide-by-zero in EPV.
## ============================================================
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
      ## Posterior expected number of coefficients whose 95%
      ## credible interval excludes zero. The fitter stores the
      ## per-coefficient quantiles in $beta$quantiles, where row
      ## 1 is the 2.5% quantile and row 5 is the 97.5% quantile
      ## (matching the user's existing BEN-Cox structure).
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


## ============================================================
## compute_epv()
##
## Events-per-variable ratio, Equation (5) of the paper:
##
##   EPV = ( sum_i delta_i^(c) ) / p_eff
##
## where delta_i^(c) is the AUGMENTED event indicator from
## Equation (4) and p_eff is the model-specific effective
## dimension from compute_p_eff().
##
## The function takes the augmented event vector directly rather
## than recomputing it. The main loop will pass in the same
## event vector that was used to fit the model, so EPV always
## reflects the actual training fold.
##
## Arguments:
##   event_train_aug  integer vector in {0,1}, the augmented
##                    event indicator delta_i^(c) of the
##                    TRAINING fold used to fit the model
##   fit              fitted-model object
##   p_raw            post-preprocessing predictor count
##
## Returns a list with:
##   $epv          numeric scalar, model-specific EPV
##   $epv_raw      numeric scalar, raw events-to-p ratio (using
##                 p_raw, not p_eff) — reported alongside per
##                 Section IV-E for transparency
##   $n_events     integer, number of events in the training fold
##   $p_eff        integer, effective dimension used in numerator
##   $p_raw        integer, raw dimension used in $epv_raw
## ============================================================
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


## ============================================================
## ph_test_global()
##
## Global proportional-hazards test for a training fold,
## implementing the diagnostic added to Section III-C in the
## revised methods section.
##
## The paper specifies: fit an unpenalized Cox model on the
## training fold and run the score test of Grambsch and
## Therneau (1994), reporting the global p-value. The test is
## descriptive — it is NOT used to gate inclusion of any model.
## Its role is to indicate, for each (dataset, regime, fold)
## cell, whether the relative performance of Cox-Time can
## plausibly be related to departures from PH rather than to
## model flexibility alone.
##
## Implementation notes:
##   - The unpenalized Cox model can fail to converge on
##     training folds with p comparable to or larger than the
##     event count (e.g., METABRIC at high censoring with
##     p_raw = 440). In such cases coxph() emits warnings about
##     "Loglik converged before variable" or "X matrix deemed
##     to be singular." We catch these and return a structured
##     NA result so the pipeline does not halt.
##   - For datasets where the Cox fit succeeds but cox.zph
##     produces NaN in the test statistic (which can happen with
##     near-zero variance covariates), we likewise return NA.
##   - When p_raw is large (say > n_events / 2), fitting the
##     unpenalized Cox is statistically meaningless even when it
##     converges, and the resulting PH test is unreliable. We
##     fit on the FULL covariate set anyway because (a) the test
##     is for the dataset, not for any specific model, and (b)
##     the alternative — pre-screening covariates — would
##     contaminate the diagnostic. The user should interpret the
##     p-value with care for high-dimensional folds, and we
##     record the n_events/p_raw ratio in the return object so
##     that downstream code can flag low-power cells.
##
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
## ============================================================
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
  
  ## Build the data frame for coxph(). Column names are taken
  ## from X_train if present, otherwise generated. Generated
  ## names use simple V1, V2, ... rather than X1, X2 to avoid
  ## collisions with the survival columns.
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