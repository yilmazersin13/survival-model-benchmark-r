## ============================================================
## 06_evaluation_metrics.R
##
## The four predictive performance metrics specified in
## Section III-C of the paper:
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
## followed by metric-specific arguments. The shape of `prediction`
## is fixed per metric:
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
## The model fitters in 05_model_fitters.R will expose
## predict_risk() and predict_survival() methods that return
## these objects directly. The main loop will call the right
## predictor for the right metric.
##
## Dependencies. We use:
##   - survival::concordance() for Harrell's C, which matches
##     the user's existing METABRIC code and is the canonical
##     R implementation.
##   - survival::survfit() for Kaplan-Meier estimates of the
##     censoring distribution, used by IBS and Uno's C.
##   - Base R for the calibration slope (a simple GLM).
##
## We deliberately do NOT use survAUC or pec for these. Both are
## reasonable alternatives but pec in particular has a heavy
## dependency footprint and survAUC is no longer well maintained.
## Reimplementing the four metrics from first principles in
## ~250 lines keeps the dependency surface small and lets us
## verify correctness against the paper's definitions directly.
## ============================================================

suppressPackageStartupMessages({
  library(survival)
})


## ============================================================
## harrell_c()
##
## Harrell's concordance index, Section III-C. Estimates the
## probability that, for a randomly chosen comparable pair, the
## subject with the higher predicted risk experiences the event
## first.
##
## We use survival::concordance() which is the canonical
## implementation. Note the sign convention: concordance() expects
## a covariate to enter Surv ~ I(...) such that LARGER values
## indicate LONGER survival, so a risk score (where larger means
## SHORTER survival) is passed in negated.
## ============================================================
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


## ============================================================
## uno_c()
##
## Uno's inverse-probability-of-censoring-weighted concordance
## measure (Uno et al., 2011), Section III-C. Mitigates the
## dependence of Harrell's C on the censoring distribution and
## is more appropriate in the high-censoring regimes studied in
## the paper.
##
## Implementation:
##   - Estimate the censoring survival function G(t) using a
##     Kaplan-Meier fit on the REVERSED event indicator
##     (1 - event), evaluated on the training data passed in
##     via `time_train` and `event_train`. Using training data
##     for G() rather than test data is the standard convention
##     and avoids contaminating the test metric with test-set
##     censoring structure.
##   - For each comparable pair (i, j) with subject i having an
##     observed event at time T_i < T_j, contribute the weight
##     1 / G(T_i)^2 to the denominator and the same weight times
##     the concordance indicator to the numerator.
##   - The truncation horizon `tau` defaults to the 90th
##     percentile of observed test event times to avoid
##     instability where G(t) approaches zero.
##
## This is implemented from first principles rather than via a
## package because (a) the standard R implementations
## (survAUC::UnoC, survC1::Inf.Cval) have inconsistent maintenance
## status, and (b) a 30-line implementation is auditable against
## Equation (4.2) of Uno et al., 2011.
## ============================================================
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
  ## Step function for G(t). approxfun with method = "constant"
  ## and f = 0 gives the right-continuous step function that is
  ## the convention for Kaplan-Meier estimators.
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


## ============================================================
## integrated_brier_score()
##
## Integrated Brier score (IBS), Section III-C. The Brier score
## at time t is estimated consistently in the presence of
## censoring using inverse-probability-of-censoring weights
## (Graf et al., 1999), and integrated over a clinically
## meaningful time horizon.
##
## Implementation follows Graf et al. (1999), Equation (12).
## For each evaluation time t and each test subject i:
##
##   contribution_i(t) =
##     1{T_i <= t, delta_i = 1} * (0 - S_hat(t | x_i))^2 / G(T_i)
##   + 1{T_i >  t}             * (1 - S_hat(t | x_i))^2 / G(t)
##
## where G(.) is the Kaplan-Meier estimator of the censoring
## distribution, fit on the training data. Subjects censored
## before time t with no event contribute zero (their fate is
## unknown).
##
## The IBS is the trapezoidal integral of BS(t) over the
## evaluation grid divided by the grid length.
##
## Arguments:
##   time_test     test fold times
##   event_test    test fold event indicators
##   prediction    n_test x length(eval_times) matrix of predicted
##                 survival probabilities
##   eval_times    numeric vector of evaluation time points,
##                 strictly increasing, within the support of
##                 the training time distribution
##   time_train    training fold times (for G estimation)
##   event_train   training fold event indicators
## ============================================================
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


## ============================================================
## calibration_slope()
##
## Calibration slope at a fixed horizon, Section III-C. The
## paper specifies: regress the observed event indicator on the
## log cumulative hazard implied by each model at a pre-specified
## evaluation horizon, defined consistently within each
## dataset-regime combination, and report the resulting slope.
##
## Implementation. The natural regression on the log cumulative
## hazard scale is a binomial GLM with the COMPLEMENTARY LOG-LOG
## link. Under cloglog link, the linear predictor is
## log(-log(1 - p)) where p is the event probability at the
## horizon, which is exactly the log cumulative hazard scale of
## a proportional-hazards model. The slope coefficient is then
## directly comparable to 1: a value of 1 indicates perfect
## calibration, values below 1 indicate overconfident risk
## predictions, and values above 1 indicate underconfident
## predictions.
##
## This is the calibration regression of Crowson, Atkinson &
## Therneau (2016, Stat Methods Med Res), and is the appropriate
## construction for fixed-horizon survival calibration. Mixing a
## cloglog-scale linear predictor with a logit link (or vice
## versa) produces a coefficient that has no clean interpretation
## and does not have target value 1.
##
## Steps:
##   - Convert predicted survival S_hat(tau | x) to predicted
##     event probability  p_hat(tau | x) = 1 - S_hat(tau | x).
##   - Form the cloglog-scale linear predictor
##       eta = log(-log(1 - p_hat)) = log(-log(S_hat)),
##     which is the log cumulative hazard at the horizon.
##   - Define the binary outcome at the horizon:
##       Z_i = 1 if T_i <= tau and delta_i = 1
##       Z_i = 0 if T_i >  tau
##     Subjects censored before tau (T_i < tau, delta_i = 0) are
##     dropped, since their status at tau is unknown.
##   - Fit  glm(Z ~ eta, family = binomial(link = "cloglog")) and
##     return the slope coefficient on eta.
##
## Arguments:
##   time_test, event_test  test fold survival data
##   prediction             numeric vector of predicted survival
##                          probabilities S_hat(tau | x), one per
##                          test subject
##   horizon                numeric scalar, the fixed evaluation
##                          time tau, in the same units as
##                          time_test (days)
## ============================================================
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


## ============================================================
## evaluate_all_metrics()
##
## Convenience wrapper that computes all four metrics for a
## single fitted model on a single test fold. Returns a one-row
## data frame in the canonical column order, suitable for
## rbind-ing across folds and models.
##
## Arguments:
##   time_test, event_test    test fold survival data
##   risk                     numeric risk score, length n_test
##                            (for harrell_c and uno_c)
##   surv_matrix              n_test x length(eval_times) matrix
##                            of predicted survival probabilities
##                            (for IBS)
##   surv_at_horizon          numeric vector of predicted S(tau | x)
##                            (for calibration_slope)
##   eval_times               IBS evaluation grid
##   horizon                  calibration evaluation horizon
##   time_train, event_train  training fold survival data, used
##                            to estimate G(.) for IBS and Uno's C
##
## Returns a data frame with columns:
##   harrell_c, uno_c, ibs, calib_slope
## ============================================================
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