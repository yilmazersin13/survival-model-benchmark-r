## ============================================================
## 05_model_fitters.R

##   fit_cox()        unpenalized Cox proportional hazards
##   fit_lasso_cox()  lasso-penalized Cox via glmnet
##   fit_rsf()        random survival forest via randomForestSRC
##   fit_ben_cox()    Bayesian elastic net Cox via Stan

##
## Every fitter takes the same core arguments:
##
##   X_train       n_train x p numeric matrix (preprocessed,
##                 from 02_preprocessing.R)
##   time_train    numeric vector of training-fold times
##   event_train   integer vector in {0,1} of training-fold
##                 event indicators
##   inner_cv      NULL or a list of inner CV splits from
##                 make_inner_cv() for hyperparameter tuning;
##                 when NULL, each fitter uses its own sensible
##                 default (e.g. cv.glmnet's built-in 10-fold)
##   seed          integer or NULL for reproducibility
##   ...           model-specific optional arguments
##
## Every fitter returns a list with these slots:
##
##   $model_class       "cox" | "lasso_cox" | "rsf" | "ben_cox"
##                      (dispatched on by compute_p_eff() in
##                      07_diagnostics.R)
##   $raw_fit           the underlying model object
##   $predict_risk      function(X_new) -> numeric vector
##                      length nrow(X_new), higher = higher hazard
##   $predict_survival  function(X_new, times) -> matrix
##                      nrow(X_new) x length(times), predicted
##                      survival probabilities S_hat(t | x)
##   $coef              numeric vector of linear coefficients
##                      (Cox / Lasso-Cox / BEN-Cox only; NULL for RSF)
##   $beta              for BEN-Cox: list(mean, quantiles) with
##                      quantiles a 5 x p matrix in the same
##                      layout as the user's existing code
##                      (rows: 2.5%, 25%, 50%, 75%, 97.5%)
##   $hyperparams       list of final hyperparameters used
##   $fit_time_sec      numeric, wall-clock seconds for the fit
##
## The predict_risk and predict_survival slots are closures that
## capture whatever state they need from the fit. Downstream code
## calls them through the uniform interface without caring which
## model is behind them.


suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(randomForestSRC)
  library(rstan)
})


.breslow_baseline <- function(time, event, lp) {
  ord <- order(time)
  t_ord <- time[ord]
  e_ord <- event[ord]
  lp_ord <- lp[ord]
  
  exp_lp <- exp(lp_ord)

  rev_cumsum <- rev(cumsum(rev(exp_lp)))
  
  ## Hazard increments at event times only.
  dH <- ifelse(e_ord == 1, 1 / rev_cumsum, 0)
  H0_at_t <- cumsum(dH)
  

  approxfun(t_ord, H0_at_t,
            method = "constant", f = 0,
            yleft = 0, yright = max(H0_at_t))
}


fit_cox <- function(X_train, time_train, event_train,
                    inner_cv = NULL, seed = NULL, ...) {
  t0 <- Sys.time()
  
  ## Build a data frame with safe column names so that the
  ## formula interface does not choke on special characters.
  col_names <- if (is.null(colnames(X_train))) {
    paste0("V", seq_len(ncol(X_train)))
  } else {
    make.names(colnames(X_train), unique = TRUE)
  }
  df <- data.frame(time = time_train, event = event_train,
                   X_train, check.names = FALSE)
  colnames(df)[-(1:2)] <- col_names
  
  formula_str <- sprintf("Surv(time, event) ~ %s",
                         paste(col_names, collapse = " + "))
  
  converged <- TRUE
  raw_fit <- withCallingHandlers(
    tryCatch(
      coxph(as.formula(formula_str), data = df,
            model = FALSE, x = FALSE, y = FALSE),
      error = function(e) NULL
    ),
    warning = function(w) {
      converged <<- FALSE
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(raw_fit)) {

    coef_vec <- rep(0, ncol(X_train))
    names(coef_vec) <- col_names
    
    out <- list(
      model_class = "cox",
      raw_fit = NULL,
      coef = coef_vec,
      predict_risk = function(X_new) rep(0, nrow(X_new)),
      predict_survival = function(X_new, times) {
        matrix(1, nrow = nrow(X_new), ncol = length(times))
      },
      hyperparams = list(converged = FALSE, degenerate = TRUE),
      fit_time_sec = as.numeric(Sys.time() - t0, units = "secs")
    )
    return(out)
  }
  
  ## Extract coefficients, defaulting NAs (perfect collinearity)
  ## to zero so that downstream predictions remain finite.
  coef_vec <- coef(raw_fit)
  coef_vec[is.na(coef_vec)] <- 0
  
  ## Baseline cumulative hazard from the training data.
  lp_train <- as.numeric(X_train %*% coef_vec)
  H0_fun <- .breslow_baseline(time_train, event_train, lp_train)
  
  predict_risk <- function(X_new) {
    as.numeric(X_new %*% coef_vec)
  }
  predict_survival <- function(X_new, times) {
    lp_new <- as.numeric(X_new %*% coef_vec)
    H0_t <- H0_fun(times)
    ## S(t|x) = exp(-H0(t) * exp(lp))
    outer(exp(lp_new), H0_t, function(r, h) exp(-h * r))
  }
  
  list(
    model_class      = "cox",
    raw_fit          = raw_fit,
    coef             = coef_vec,
    predict_risk     = predict_risk,
    predict_survival = predict_survival,
    hyperparams      = list(converged = converged),
    fit_time_sec     = as.numeric(Sys.time() - t0, units = "secs")
  )
}



fit_lasso_cox <- function(X_train, time_train, event_train,
                          inner_cv = NULL, seed = NULL, ...) {
  t0 <- Sys.time()
  if (!is.null(seed)) set.seed(seed)
  
  y <- Surv(time_train, event_train)
  

  foldid <- NULL
  if (!is.null(inner_cv)) {
    n_train <- nrow(X_train)

  }
  
  cv_fit <- cv.glmnet(x = X_train, y = y, family = "cox",
                      alpha = 1, foldid = foldid)
  raw_fit <- glmnet(x = X_train, y = y, family = "cox",
                    alpha = 1, lambda = cv_fit$lambda.min)
  
  coef_vec <- as.numeric(coef(raw_fit))
  names(coef_vec) <- colnames(X_train)
  
  lp_train <- as.numeric(X_train %*% coef_vec)
  H0_fun <- .breslow_baseline(time_train, event_train, lp_train)
  
  predict_risk <- function(X_new) {
    as.numeric(X_new %*% coef_vec)
  }
  predict_survival <- function(X_new, times) {
    lp_new <- as.numeric(X_new %*% coef_vec)
    H0_t <- H0_fun(times)
    outer(exp(lp_new), H0_t, function(r, h) exp(-h * r))
  }
  
  list(
    model_class      = "lasso_cox",
    raw_fit          = raw_fit,
    coef             = coef_vec,
    predict_risk     = predict_risk,
    predict_survival = predict_survival,
    hyperparams = list(
      lambda_min = cv_fit$lambda.min,
      lambda_1se = cv_fit$lambda.1se,
      n_nonzero  = sum(coef_vec != 0)
    ),
    fit_time_sec = as.numeric(Sys.time() - t0, units = "secs")
  )
}



fit_rsf <- function(X_train, time_train, event_train,
                    inner_cv = NULL, seed = NULL,
                    ntree = 500, ...) {
  t0 <- Sys.time()
  if (!is.null(seed)) set.seed(seed)
  
  df <- data.frame(time = time_train, event = event_train, X_train,
                   check.names = TRUE)
  

  tuned <- tryCatch(
    tune.rfsrc(Surv(time, event) ~ ., data = df,
               mtryStart = max(1, floor(sqrt(ncol(X_train)))),
               nodesizeTry = c(5, 10, 15, 25),
               ntreeTry = 100,
               trace = FALSE),
    error = function(e) NULL
  )
  if (is.null(tuned)) {
    best_mtry <- max(1, floor(sqrt(ncol(X_train))))
    best_nodesize <- 15
  } else {
    best_mtry     <- tuned$optimal["mtry"]
    best_nodesize <- tuned$optimal["nodesize"]
  }
  
  raw_fit <- rfsrc(Surv(time, event) ~ ., data = df,
                   ntree = ntree,
                   mtry = best_mtry,
                   nodesize = best_nodesize,
                   importance = FALSE,
                   forest = TRUE)
  

  base_times <- raw_fit$time.interest
  
  predict_risk <- function(X_new) {
    df_new <- data.frame(X_new, check.names = TRUE)
    pred <- predict.rfsrc(raw_fit, newdata = df_new,
                          importance = "none")

    as.numeric(pred$predicted)
  }
  
  predict_survival <- function(X_new, times) {
    df_new <- data.frame(X_new, check.names = TRUE)
    pred <- predict.rfsrc(raw_fit, newdata = df_new,
                          importance = "none")

    base_surv <- pred$survival
    n_new <- nrow(base_surv)
    out <- matrix(1, nrow = n_new, ncol = length(times))
    for (k in seq_along(times)) {
      ## Index of the largest base_time <= times[k]
      idx <- findInterval(times[k], base_times)
      if (idx == 0) {
        out[, k] <- 1
      } else {
        out[, k] <- base_surv[, idx]
      }
    }
    out
  }
  
  list(
    model_class      = "rsf",
    raw_fit          = raw_fit,
    coef             = NULL,
    predict_risk     = predict_risk,
    predict_survival = predict_survival,
    hyperparams = list(
      ntree    = ntree,
      mtry     = as.integer(best_mtry),
      nodesize = as.integer(best_nodesize)
    ),
    fit_time_sec = as.numeric(Sys.time() - t0, units = "secs")
  )
}


.ben_cox_stan_code <- "
data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  vector[N] y;
  int<lower=0, upper=1> event[N];
}
parameters {
  vector[P] beta_raw;
  vector<lower=1e-6>[P] tau_sq;
  real<lower=0.01> lambda1;
  real<lower=0.01> lambda2;
}
transformed parameters {
  vector[P] beta;
  for (j in 1:P) {
    real variance_j = tau_sq[j] / (1.0 + lambda2 * tau_sq[j]);
    beta[j] = beta_raw[j] * sqrt(variance_j);
  }
}
model {
  lambda1 ~ cauchy(0, 1);
  lambda2 ~ cauchy(0, 1);
  for (j in 1:P) {
    tau_sq[j] ~ exponential(0.5 * square(lambda1));
  }
  beta_raw ~ std_normal();
  {
    vector[N] risk = X * beta;
    real current_log_sum = negative_infinity();
    real log_lik = 0;
    for (i in 1:N) {
      int idx = N - i + 1;
      current_log_sum = log_sum_exp(current_log_sum, risk[idx]);
      if (event[idx] == 1) {
        log_lik += risk[idx] - current_log_sum;
      }
    }
    target += log_lik;
  }
}
"



fit_ben_cox <- function(X_train, time_train, event_train,
                        inner_cv = NULL, seed = NULL,
                        n_iter = 1500, n_warmup = 750,
                        n_chains = 2, adapt_delta = 0.90, ...) {
  t0 <- Sys.time()
  if (!is.null(seed)) set.seed(seed)
  
  ## Stan likelihood uses a backward cumulative sum that assumes
  ## observations are sorted by time in ascending order.
  ord <- order(time_train, decreasing = FALSE)
  X_sorted <- X_train[ord, , drop = FALSE]
  t_sorted <- time_train[ord]
  e_sorted <- as.integer(event_train[ord])
  
  stan_data <- list(
    N = nrow(X_sorted),
    P = ncol(X_sorted),
    X = X_sorted,
    y = t_sorted,
    event = e_sorted
  )
  
  raw_fit <- rstan::stan(
    model_code = .ben_cox_stan_code,
    data       = stan_data,
    iter       = n_iter,
    warmup     = n_warmup,
    chains     = n_chains,
    control    = list(adapt_delta = adapt_delta, max_treedepth = 12),
    seed       = if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed,
    refresh    = 0,
    verbose    = FALSE,
    init = function() list(
      beta_raw = rnorm(ncol(X_sorted), 0, 0.1),
      tau_sq   = rep(0.1, ncol(X_sorted)),
      lambda1  = 1,
      lambda2  = 1
    )
  )
  
  samples <- rstan::extract(raw_fit)
  beta_mean <- colMeans(samples$beta)
  ## Quantile matrix in the 5 x P layout expected by
  ## compute_p_eff(): rows 2.5%, 25%, 50%, 75%, 97.5%.
  beta_q <- apply(samples$beta, 2, quantile,
                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  lp_train <- as.numeric(X_train %*% beta_mean)
  H0_fun <- .breslow_baseline(time_train, event_train, lp_train)
  
  predict_risk <- function(X_new) {
    as.numeric(X_new %*% beta_mean)
  }
  predict_survival <- function(X_new, times) {
    lp_new <- as.numeric(X_new %*% beta_mean)
    H0_t <- H0_fun(times)
    outer(exp(lp_new), H0_t, function(r, h) exp(-h * r))
  }
  
  list(
    model_class      = "ben_cox",
    raw_fit          = raw_fit,
    coef             = beta_mean,
    beta             = list(mean = beta_mean, quantiles = beta_q),
    predict_risk     = predict_risk,
    predict_survival = predict_survival,
    hyperparams = list(
      n_iter      = n_iter,
      n_warmup    = n_warmup,
      n_chains    = n_chains,
      adapt_delta = adapt_delta,
      lambda1_mean = mean(samples$lambda1),
      lambda2_mean = mean(samples$lambda2)
    ),
    fit_time_sec = as.numeric(Sys.time() - t0, units = "secs")
  )
}
