## ============================================================
## 03_censoring_augmentation.R

## Internal helper: given a candidate threshold tau, draw the
## auxiliary censoring times and return the resulting augmented
## (time, event) and the empirical censoring rate.
##
## This is factored out of the bisection loop so that the same
## logic is used both inside the search and for the final draw
## with the chosen tau.

.augment_with_tau <- function(time, event, tau, rng_state = NULL) {
  if (!is.null(rng_state)) {

    assign(".Random.seed", rng_state, envir = .GlobalEnv)
  }
  n <- length(time)
  C_star <- runif(n, min = 0, max = tau)
  
  time_aug  <- pmin(time, C_star)
  event_aug <- as.integer(event * (time <= C_star))
  
  cens_rate <- 1 - mean(event_aug)
  list(time = time_aug, event = event_aug, cens_rate = cens_rate)
}



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
  

  t_max <- max(time)
  tau_lo <- t_max * 1e-4
  tau_hi <- t_max * 100
  

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
