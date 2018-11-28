# some helper functions we'll use throughout

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

# compute log of raw importance ratios
# sums over observations *not* over posterior samples
sum_log_ratios <- function(ll, ids = NULL) {
  if (!is.null(ids)) ll <- ll[, ids, drop = FALSE]
  - rowSums(ll)
}

# for printing comparisons later
rbind_print <- function(...) {
  round(rbind(...), digits = 2)
}

scale_unit_interval <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

summarize_elpds <- function(elpds) {
  elpds <- na.omit(elpds)
  c(ELPD = sum(elpds), SE = sqrt(length(elpds) * var(elpds)))
}

seq_pos <- function(from, to, ...) {
  if (from > to) {
    out <- integer(0)
  } else {
    out <- seq.int(from, to, ...)
  }
  out
}

# Perform exact leave-future-out cross-validation
# 
# @param fit a brmsfit object
# @param M steps to predict into the future
# @param L minimal number of observations required for fitting the model
# 
# @return A vector of pointwise ELPD values
exact_lfo <- function(fit, M, L, B = NA) {
  require(brms)
  stopifnot(is.brmsfit(fit))
  df <- model.frame(fit)
  N <- NROW(df)
  if (!length(B) || anyNA(B)) {
    B <- N
  }
  if (B < M) {
    stop2("The left-out block must at least consist ", 
          "of all predicted observations")
  }
  
  # compute exact LFO likelihoods
  ids <- (N - M + 1):max(L + 1, 2)
  loglikm <- matrix(nrow = nsamples(fit), ncol = N)
  for (i in ids) {
    ioos <- 1:(i + M - 1)
    oos <- i:(i + M - 1)
    ind_rm <- i:min(i + B - 1, N)
    newdf <- df
    # future and past values should not be treated as one time-series
    # in the block version as they have a block of missing obs between them
    newdf$group <- c(rep("past", i - 1), rep("future", nrow(newdf) - i + 1))
    newdf <- newdf[-ind_rm, , drop = FALSE]
    fit_i <- update(fit, newdata = newdf, recompile = FALSE, refresh = 0)
    ll <- log_lik(fit_i, newdata = df[ioos, , drop = FALSE], oos = oos)
    loglikm[, i] <- rowSums(ll[, oos, drop = FALSE])
  }
  # compute and return elpds per observation
  apply(loglikm, 2, log_mean_exp)
}

approx_lfo <- function(fit, M, L, B = NA, k_thres = 0.6) {
  require(brms)
  require(loo)
  stopifnot(is.brmsfit(fit))
  df <- model.frame(fit)
  N <- NROW(df)
  if (!length(B) || anyNA(B)) {
    B <- Inf
  }
  if (B < M) {
    stop2("The left-out block must at least consist ", 
          "of all predicted observations (B >= M).")
  }
  
  # compute approximate LFO likelihoods
  loglikm <- loglik <- matrix(nrow = nsamples(fit), ncol = N)
  out <- rep(NA, N)
  fit_i <- fit
  # observations to predict
  ids <- (N - M + 1):(L + 1)
  # last observation included in the model fitting
  i_refit <- N
  refits <- ks <- numeric(0)
  # no isolated predictions of the last M - 1 observations
  ind_init <- seq_pos(N - M + 2, N)
  if (length(ind_init)) {
    loglik[, ind_init] <- log_lik(fit_i)[, ind_init, drop = FALSE] 
  }
  for (i in ids) {
    ioos <- 1:(i + M - 1)
    oos <- i:(i + M - 1)
    ll <- log_lik(fit_i, newdata = df[ioos, , drop = FALSE], oos = oos)
    loglikm[, i] <- rowSums(ll[, oos, drop = FALSE])
    loglik[, i] <- ll[, i]
    # observations over which to perform importance sampling
    ind1_logratio <- i:min(i + B - 1, i_refit)
    logratio <- sum_log_ratios(loglik, ind1_logratio)
    if (B < Inf) {
      # in the block version some observations need to be added back again
      ind2_logratio <- seq_pos(max(i + B, i_refit + 1), min(i_refit + B, N))
      logratio <- logratio - sum_log_ratios(loglik, ind2_logratio)
    }
    psis_part <- suppressWarnings(psis(logratio))
    k <- pareto_k_values(psis_part)
    ks <- c(ks, k)
    if (k > k_thres) {
      # refit the model based on the first i - 1 observations
      i_refit <- i - 1
      refits <- c(refits, i)
      ind_rm <- i:min(i + B - 1, N)
      newdf <- df
      # future and past values should not be treated as one time-series
      # in the block version as they have a block of missing obs between them
      newdf$group <- c(rep("past", i - 1), rep("future", nrow(newdf) - i + 1))
      newdf <- newdf[-ind_rm, , drop = FALSE]
      fit_i <- update(fit, newdata = newdf, recompile = FALSE, refresh = 0)
      # perform exact LFO for the ith observation
      ll <- log_lik(fit_i, newdata = df[ioos, , drop = FALSE], oos = oos)
      loglik[, i] <- ll[, i]
      loglikm[, i] <- rowSums(ll[, oos, drop = FALSE])
      out[i] <- log_mean_exp(loglikm[, i])
    } else {
      # PSIS approximate LFO is possible
      lw_i <- weights(psis_part, normalize = TRUE)[, 1]
      out[i] <- log_sum_exp(lw_i + loglikm[, i])
    }
  }
  attr(out, "ks") <- ks
  attr(out, "refits") <- refits
  out
}

fit_model <- function(cond, ...) {
  require(brms)
  stopifnot(is.data.frame(cond) && NROW(cond) == 1L)
  N <- cond$N
  model <- cond$model
  time <- seq_len(N)
  stime <- scale_unit_interval(time)
  df <- data.frame(time = time, stime = stime, group = "past")
  ar_prior <- prior(normal(0, 0.5), class = "ar")
  ar_autocor <- cor_ar(~ stime | group, p = 2)
  if (model == "constant") {
    df$y <- rnorm(N)
    fit <- brm(y ~ 1, data = df, refresh = 0, ...)
  } else if (model == "linear") {
    df$y <- 17 * stime + rnorm(N) 
    fit <- brm(y ~ stime, data = df, refresh = 0, ...)
  } else if (model == "quadratic") {
    df$y <- 17 * stime - 25 * stime^2 + rnorm(N)
    fit <- brm(y ~ stime + I(stime^2), data = df, refresh = 0, ...)
  } else if (model == "AR2_only") {
    df$y <- as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y ~ 1, data = df, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else if (model == "AR2_linear") {
    df$y <- 17 * stime +
      as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y ~ stime, data = df, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else if (model == "AR2_quadratic") {
    df$y <- 17 * stime - 25 * stime^2 + 
      as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y ~ stime + I(stime^2), 
      data = df, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else {
    stop("Model '", model, "' is not supported.")
  }
  fit
}

sim_fun <- function(j, conditions, ...) {
  # simulate data and fit the corresponding model
  fit <- fit_model(cond = conditions[j, ], ...)
  N <- conditions$N[j]
  M <- conditions$M[j]
  L <- conditions$L[j]
  B <- conditions$B[j]
  k_thres <- conditions$k_thres[j]
  # LOO-CV
  loo_cv <- loo(log_lik(fit)[, (L + 1):N])
  # Exact LFO-CV
  lfo_exact_elpds <- exact_lfo(fit, M = M, L = L, B = B)
  lfo_exact_elpd <- summarize_elpds(lfo_exact_elpds)
  # Approximate LFO-CV
  lfo_approx_elpds <- approx_lfo(
    fit, M = M, L = L, B = B, k_thres = k_thres
  )
  lfo_approx_elpd <- summarize_elpds(lfo_approx_elpds)
  # return all relevant information
  list(
    # fit requires too much space
    loo_cv = loo_cv,
    lfo_exact_elpds = lfo_exact_elpds,
    lfo_exact_elpd = lfo_exact_elpd,
    lfo_approx_elpds = lfo_approx_elpds,
    lfo_approx_elpd = lfo_approx_elpd
  )
}
