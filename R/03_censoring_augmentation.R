## ============================================================
## 03_censoring_augmentation.R
##
## Censoring augmentation protocol from Section IV-C of the paper.
##
## Implements Equation (4):
##
##   T_tilde_i^(c) = min(T_tilde_i, C_i_star)
##   delta_i^(c)   = delta_i * 1{T_tilde_i <= C_i_star}
##
## where C_i_star ~ Uniform(0, tau_c) independently of T_i and
## x_i, and tau_c is chosen by bisection so that the empirical
## censoring rate of the augmented data matches a target c up to
## a tolerance of 0.01.
##
## Two important properties from Section IV-C, both enforced here:
##
##   1. The procedure can only INCREASE censoring relative to the
##      baseline. If the requested target c is less than or equal
##      to the baseline censoring rate, augmentation is a no-op
##      and the function returns the input unchanged with a
##      message. This is the case for TCGA-BRCA at any c <= 0.85
##      under the present datasets.
##
##   2. Non-informative censoring is preserved by construction
##      because C_i_star is drawn independently of T_i and x_i.
##      The covariate matrix X is therefore not touched.
##
## The function operates on a single fold (typically the training
## fold of an outer CV split) and is called from the main
## benchmark loop in 04_cv_scaffolding.R / main.R, NOT at load
## time. This matches the design of Section IV-C, where
## augmentation is applied per training fold.
## ============================================================

## ------------------------------------------------------------
## Internal helper: given a candidate threshold tau, draw the
## auxiliary censoring times and return the resulting augmented
## (time, event) and the empirical censoring rate.
##
## This is factored out of the bisection loop so that the same
## logic is used both inside the search and for the final draw
## with the chosen tau.
## ------------------------------------------------------------
.augment_with_tau <- function(time, event, tau, rng_state = NULL) {
  if (!is.null(rng_state)) {
    ## Restore the RNG state so that, given a fixed seed and a
    ## fixed tau, the augmented data are bit-for-bit reproducible.
    assign(".Random.seed", rng_state, envir = .GlobalEnv)
  }
  n <- length(time)
  C_star <- runif(n, min = 0, max = tau)
  
  time_aug  <- pmin(time, C_star)
  event_aug <- as.integer(event * (time <= C_star))
  
  cens_rate <- 1 - mean(event_aug)
  list(time = time_aug, event = event_aug, cens_rate = cens_rate)
}


## ============================================================
## augment_censoring()
##
## Apply the censoring augmentation protocol of Section IV-C to
## a single (time, event) pair, targeting a specified censoring
## rate `c_target`.
##
## Arguments:
##   time       numeric vector, T_tilde_i (the observed time,
##              already min(T_i, C_i) at the dataset level)
##   event      integer vector in {0,1}, delta_i
##   c_target   numeric scalar in (0,1), the desired censoring
##              proportion after augmentation
##   tolerance  numeric scalar, max allowed |achieved - target|;
##              default 0.01 to match Section IV-C
##   max_iter   integer, max bisection iterations; default 60
##              (more than enough to hit 0.01 tolerance)
##   seed       integer or NULL; if non-NULL, sets the RNG seed
##              before drawing C_star, so that the augmented
##              data are reproducible across runs
##
## Returns a list with:
##   $time          augmented T_tilde_i^(c), length n
##   $event         augmented delta_i^(c), length n
##   $c_target      the requested target
##   $c_achieved    the empirical censoring rate after augmentation
##   $tau           the threshold tau_c selected by bisection
##   $augmented     logical; FALSE if c_target <= baseline (no-op)
## ============================================================
augment_censoring <- function(time, event,
                              c_target,
                              tolerance = 0.01,
                              max_iter = 60,
                              seed = NULL) {
  ## --- Argument checks ---
  if (length(time) != length(event)) {
    stop("time and event must have the same length")
  }
  if (!all(event %in% c(0L, 1L, 0, 1))) {
    stop("event must be 0/1")
  }
  if (any(time <= 0)) {
    stop("time must be strictly positive")
  }
  if (c_target <= 0 || c_target >= 1) {
    stop("c_target must be in the open interval (0, 1)")
  }
  
  baseline_cens <- 1 - mean(event)
  
  ## --- Property 1 from Section IV-C: only increase censoring ---
  ## If the target is at or below the baseline, return the input
  ## unchanged. The paper handles this case explicitly for TCGA-BRCA.
  if (c_target <= baseline_cens + tolerance) {
    message(sprintf(
      "augment_censoring: target c = %.3f <= baseline %.3f; ",
      c_target, baseline_cens),
      "returning input unchanged (no augmentation applied).")
    return(list(
      time       = time,
      event      = as.integer(event),
      c_target   = c_target,
      c_achieved = baseline_cens,
      tau        = NA_real_,
      augmented  = FALSE
    ))
  }
  
  ## --- Bisection over tau ---
  ## Monotonicity: as tau decreases, more events get cut off by
  ## the auxiliary censoring time, so the censoring rate increases.
  ## Therefore tau is searched on the interval (tau_lo, tau_hi)
  ## with the convention that small tau -> high censoring.
  ##
  ## Lower bound: a very small fraction of the maximum observed
  ## time. Upper bound: well beyond the maximum, so that the
  ## auxiliary censoring is essentially never the minimum and
  ## the augmented censoring rate is close to baseline.
  t_max <- max(time)
  tau_lo <- t_max * 1e-4
  tau_hi <- t_max * 100
  
  ## Set the seed once before the search so that every candidate
  ## tau is evaluated against the *same* random draws of C_star.
  ## Without this, the search target moves under the bisection
  ## and convergence becomes erratic. This is a standard trick
  ## for stochastic root-finding.
  if (!is.null(seed)) set.seed(seed)
  rng_state <- get(".Random.seed", envir = .GlobalEnv)
  
  achieved <- NA_real_
  tau_mid  <- NA_real_
  
  for (iter in seq_len(max_iter)) {
    tau_mid <- 0.5 * (tau_lo + tau_hi)
    aug <- .augment_with_tau(time, event, tau_mid, rng_state = rng_state)
    achieved <- aug$cens_rate
    
    if (abs(achieved - c_target) <= tolerance) {
      break
    }
    if (achieved > c_target) {
      ## Too much censoring -> tau too small -> raise tau
      tau_lo <- tau_mid
    } else {
      ## Too little censoring -> tau too large -> lower tau
      tau_hi <- tau_mid
    }
  }
  
  if (abs(achieved - c_target) > tolerance) {
    warning(sprintf(
      "augment_censoring: bisection did not converge within %d iterations. ",
      max_iter),
      sprintf("Target = %.3f, achieved = %.3f, tau = %.3g.",
              c_target, achieved, tau_mid),
      call. = FALSE)
  }
  
  ## --- Final augmentation with the chosen tau ---
  ## Using the same rng_state guarantees the returned (time, event)
  ## are exactly the data that achieved the reported censoring rate.
  final <- .augment_with_tau(time, event, tau_mid, rng_state = rng_state)
  
  list(
    time       = final$time,
    event      = final$event,
    c_target   = c_target,
    c_achieved = final$cens_rate,
    tau        = tau_mid,
    augmented  = TRUE
  )
}


## ============================================================
## augment_dataset()
##
## Convenience wrapper that takes a full dataset object (the
## kind returned by the loaders in 01_data_loaders.R) and
## returns a new dataset object with the time and event slots
## replaced by their augmented counterparts. The X matrix and
## all metadata slots are passed through unchanged, since
## Section IV-C augments only the censoring structure.
##
## A new slot `$augmentation` is added recording the target,
## achieved rate, tau, and whether augmentation was actually
## applied. This is what the EPV computation in §IV-E and the
## results tables in §V will read from.
##
## Arguments:
##   dataset    a dataset object from 01_data_loaders.R
##   c_target   target censoring rate
##   ...        passed through to augment_censoring()
## ============================================================
augment_dataset <- function(dataset, c_target, ...) {
  required <- c("time", "event", "X", "dataset_name", "baseline_censoring")
  if (!all(required %in% names(dataset))) {
    stop("augment_dataset: input does not look like a loader object")
  }
  
  aug <- augment_censoring(dataset$time, dataset$event,
                           c_target = c_target, ...)
  
  out <- dataset
  out$time  <- aug$time
  out$event <- aug$event
  out$augmentation <- list(
    c_target   = aug$c_target,
    c_achieved = aug$c_achieved,
    tau        = aug$tau,
    augmented  = aug$augmented
  )
  out
}