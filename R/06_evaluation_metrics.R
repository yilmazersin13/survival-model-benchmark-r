## 06_evaluation_metrics.R
##

##
##   1. Harrell's concordance index               -> harrell_c()
##   2. Uno's time-dependent concordance          -> uno_c()
##   3. Integrated Brier score                    -> integrated_brier_score()
##   4. Calibration slope at a fixed horizon      -> calibration_slope()
##
## Calling convention. Every metric takes the same first three
## arguments:
##
##   time_test   numeric vector, observed times in the test fold
##   event_test  integer vector in {0,1}, event indicators
##   prediction  the prediction object (shape depends on metric)
##

##
##   harrell_c              numeric risk score, one per subject
##                          (higher = higher predicted hazard)
##   uno_c                  numeric risk score, one per subject
##   integrated_brier_score numeric matrix, n_test x length(times),
##                          containing predicted survival
##                          probabilities S(t | x) at each
##                          evaluation time
##   calibration_slope      numeric vector of length n_test,
##                          predicted survival probability at the
##                          fixed horizon, one per subject
##

## Dependencies. We use:
##   - survival::concordance() for Harrell's C, which matches
##     the user's existing METABRIC code and is the canonical
##     R implementation.
##   - survival::survfit() for Kaplan-Meier estimates of the
##     censoring distribution, used by IBS and Uno's C.
##   - Base R for the calibration slope (a simple GLM).



suppressPackageStartupMessages({
  library(survival)
})



harrell_c <- function(time_test, event_test, prediction) {
  if (length(prediction) != length(time_test)) {
    stop("harrell_c: prediction length does not match time_test")
  }
  if (any(!is.finite(prediction))) {
    warning("harrell_c: non-finite values in prediction, returning NA",
            call. = FALSE)
    return(NA_real_)
  }
  surv_df <- data.frame(time = time_test, event = event_test)
  fit <- tryCatch(
    survival::concordance(survival::Surv(time, event) ~ I(-prediction),
                          data = surv_df),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  unname(fit$concordance)
}



uno_c <- function(time_test, event_test, prediction,
                  time_train, event_train, tau = NULL) {
  if (length(prediction) != length(time_test)) {
    stop("uno_c: prediction length does not match time_test")
  }
  if (missing(time_train) || missing(event_train)) {
    stop("uno_c: time_train and event_train are required ",
         "to estimate the censoring distribution")
  }
  if (any(!is.finite(prediction))) {
    warning("uno_c: non-finite values in prediction, returning NA",
            call. = FALSE)
    return(NA_real_)
  }
  
  ## --- Truncation horizon ---
  if (is.null(tau)) {
    event_times_test <- time_test[event_test == 1]
    if (length(event_times_test) < 5) return(NA_real_)
    tau <- quantile(event_times_test, 0.90)
  }
  
  ## --- Censoring KM on training data (reversed indicator) ---
  km_cens <- survival::survfit(
    survival::Surv(time_train, 1 - event_train) ~ 1
  )

  G_fun <- approxfun(km_cens$time, km_cens$surv,
                     method = "constant", f = 0,
                     yleft = 1, yright = min(km_cens$surv))
  
  n <- length(time_test)
  num <- 0
  den <- 0
  for (i in seq_len(n)) {
    if (event_test[i] != 1) next
    if (time_test[i] >= tau) next
    Gi <- G_fun(time_test[i])
    if (!is.finite(Gi) || Gi <= 0) next
    w <- 1 / (Gi * Gi)
    for (j in seq_len(n)) {
      if (j == i) next
      if (time_test[j] <= time_test[i]) next   # j must outlive i
      den <- den + w
      if (prediction[i] > prediction[j]) {
        num <- num + w
      } else if (prediction[i] == prediction[j]) {
        num <- num + 0.5 * w
      }
    }
  }
  if (den == 0) return(NA_real_)
  num / den
}



integrated_brier_score <- function(time_test, event_test, prediction,
                                   eval_times,
                                   time_train, event_train) {
  if (!is.matrix(prediction)) {
    stop("integrated_brier_score: prediction must be a matrix")
  }
  if (nrow(prediction) != length(time_test)) {
    stop("integrated_brier_score: prediction has wrong number of rows")
  }
  if (ncol(prediction) != length(eval_times)) {
    stop("integrated_brier_score: prediction has wrong number of columns")
  }
  if (any(diff(eval_times) <= 0)) {
    stop("integrated_brier_score: eval_times must be strictly increasing")
  }
  if (any(!is.finite(prediction))) {
    warning("integrated_brier_score: non-finite predictions, returning NA",
            call. = FALSE)
    return(NA_real_)
  }
  
  ## --- Censoring KM on training data ---
  km_cens <- survival::survfit(
    survival::Surv(time_train, 1 - event_train) ~ 1
  )
  G_fun <- approxfun(km_cens$time, km_cens$surv,
                     method = "constant", f = 0,
                     yleft = 1, yright = min(km_cens$surv))
  
  n_test <- length(time_test)
  n_t    <- length(eval_times)
  
  bs <- numeric(n_t)
  for (k in seq_len(n_t)) {
    t_k <- eval_times[k]
    G_t <- G_fun(t_k)
    if (!is.finite(G_t) || G_t <= 0) {
      bs[k] <- NA_real_
      next
    }
    s_hat <- prediction[, k]
    contrib <- numeric(n_test)
    for (i in seq_len(n_test)) {
      if (time_test[i] <= t_k && event_test[i] == 1) {
        Gi <- G_fun(time_test[i])
        if (!is.finite(Gi) || Gi <= 0) next
        contrib[i] <- (0 - s_hat[i])^2 / Gi
      } else if (time_test[i] > t_k) {
        contrib[i] <- (1 - s_hat[i])^2 / G_t
      }
      ## else: censored before t_k with no event, contributes 0
    }
    bs[k] <- mean(contrib)
  }
  
  if (all(is.na(bs))) return(NA_real_)
  
  ## Trapezoidal integration over eval_times, normalized by
  ## the time-grid length so that IBS is comparable across
  ## different choices of grid.
  valid <- !is.na(bs)
  if (sum(valid) < 2) return(NA_real_)
  bs_v <- bs[valid]
  t_v  <- eval_times[valid]
  area <- sum(diff(t_v) * (bs_v[-1] + bs_v[-length(bs_v)]) / 2)
  area / (max(t_v) - min(t_v))
}



calibration_slope <- function(time_test, event_test, prediction,
                              horizon) {
  if (length(prediction) != length(time_test)) {
    stop("calibration_slope: prediction length does not match time_test")
  }
  if (any(!is.finite(prediction))) {
    warning("calibration_slope: non-finite predictions, returning NA",
            call. = FALSE)
    return(NA_real_)
  }
  if (horizon <= 0) stop("calibration_slope: horizon must be positive")
  
  ## Drop subjects censored before the horizon: their status at
  ## tau is unknown so they cannot enter the calibration target.
  keep <- !(time_test < horizon & event_test == 0)
  t_k <- time_test[keep]
  e_k <- event_test[keep]
  s_k <- prediction[keep]
  
  if (length(t_k) < 10) return(NA_real_)
  
  ## Bound away from 0 and 1 so that log(-log(.)) is finite.
  s_k <- pmin(pmax(s_k, 1e-6), 1 - 1e-6)
  
  ## Log cumulative hazard at the horizon. With a cloglog-link
  ## binomial GLM, this is the natural linear-predictor scale.
  eta <- log(-log(s_k))
  
  ## Binary outcome at the horizon.
  Z <- as.integer(t_k <= horizon & e_k == 1)
  
  if (length(unique(Z)) < 2) return(NA_real_)
  if (var(eta) < .Machine$double.eps) return(NA_real_)
  
  fit <- tryCatch(
    glm(Z ~ eta, family = binomial(link = "cloglog")),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(fit)) return(NA_real_)
  unname(coef(fit)["eta"])
}


evaluate_all_metrics <- function(time_test, event_test,
                                 risk,
                                 surv_matrix,
                                 surv_at_horizon,
                                 eval_times,
                                 horizon,
                                 time_train, event_train) {
  data.frame(
    harrell_c = harrell_c(time_test, event_test, risk),
    uno_c     = uno_c(time_test, event_test, risk,
                      time_train = time_train,
                      event_train = event_train),
    ibs       = integrated_brier_score(time_test, event_test,
                                       surv_matrix,
                                       eval_times = eval_times,
                                       time_train = time_train,
                                       event_train = event_train),
    calib_slope = calibration_slope(time_test, event_test,
                                    surv_at_horizon,
                                    horizon = horizon)
  )
}
