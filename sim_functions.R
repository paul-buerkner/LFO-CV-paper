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
  c(Estimate = sum(elpds), SE = sqrt(length(elpds) * var(elpds)))
}

plot_ks <- function(ks, ids, thres = 0.6) {
  dat_ks <- data.frame(ks = ks, ids = ids)
  ggplot(dat_ks, aes(x = ids, y = ks)) + 
    geom_point(aes(color = ks > thres), shape = 3, show.legend = FALSE) + 
    geom_hline(yintercept = thres, linetype = 2, color = "red2") + 
    scale_color_manual(values = c("cornflowerblue", "darkblue")) + 
    labs(x = "Data point", y = "Pareto k") + 
    ylim(-0.5, 1.5)
}

seq_pos <- function(from, to, ...) {
  if (from > to) {
    out <- integer(0)
  } else {
    out <- seq.int(from, to, ...)
  }
  out
}

rmse <- function(x, y, weights = NULL) {
  # args:
  #   x: SN matrix of posterior predictions
  #   y: N vector of reponse values
  #   weights: optional S vector of weights
  x <- as.matrix(x)
  stopifnot(length(y) == NCOL(x))
  out <- 0
  if (is.null(weights)) {
    for (i in seq_along(y)) {
      out <- out + sqrt(mean((x[, i] - y[i])^2))
    }
  } else {
    stopifnot(length(weights) == NROW(x))
    for (i in seq_along(y)) {
      out <- out + sqrt(sum(weights * (x[, i] - y[i])^2))
    }
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
exact_lfo <- function(fit, M, L, B = NA, criterion = c("elpd", "rmse")) {
  require(brms)
  stopifnot(is.brmsfit(fit))
  criterion <- match.arg(criterion)
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
  out <- rep(NA, N)
  conv <- vector("list", N)
  loglikm <- matrix(nrow = nsamples(fit), ncol = N)
  for (i in (N - M):L) {
    ioos <- 1:(i + M)
    oos <- (i + 1):(i + M)
    ind_rm <- (i + 1):min(i + B, N)
    newdf <- df
    if (max(ind_rm) < N) {
      # responses of the left-out block need to be coded as 
      # missing to retain a single consistent time series
      newdf$y[ind_rm] <- NA
    } else {
      newdf <- newdf[-ind_rm, , drop = FALSE]
    }
    fit_star <- update(fit, newdata = newdf, recompile = FALSE, refresh = 0)
    conv[[i]] <- convergence_summary(fit_star)
    out[i] <- lfo_criterion(
      fit_star, data = df[ioos, , drop = FALSE], oos = oos, 
      criterion = criterion
    )
  }
  attr(out, "conv") <- conv
  out
}

approx_lfo <- function(fit, M, L, B = NA, k_thres = 0.6, 
                       criterion = c("elpd", "rmse")) {
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
  out <- ks <- rep(NA, N)
  conv <- vector("list", N)
  fit_star <- fit
  # last observation included in the model fitting
  i_star <- N
  refits <- numeric(0)
  # no isolated predictions of the last M observations
  if (M > 1) {
    loglik[, (N - M + 2):N] <- log_lik(fit_star)[, (N - M + 2):N, drop = FALSE]  
  }
  for (i in (N - M):L) {
    ioos <- 1:(i + M)
    oos <- (i + 1):(i + M)
    ll <- log_lik(fit_star, newdata = df[ioos, , drop = FALSE], oos = oos)
    loglik[, i + 1] <- ll[, i + 1]
    # observations over which to perform importance sampling
    ilr1 <- (i + 1):min(i + B, i_star)
    logratio <- sum_log_ratios(loglik, ilr1)
    if (B < Inf) {
      # in the block version some observations need to be added back again
      ilr2 <- seq_pos(max(i + B + 1, i_star + 1), min(i_star + B, N))
      if (length(ilr2)) {
        # observations in the left-out block are modeled as missing
        ind_B <- (i + 1):(i + B)
        subdf <- df[seq_len(max(ilr2)), , drop = FALSE]
        ll_after_block <- log_lik(fit_star, newdata = subdf, oos = ind_B)
        logratio <- logratio - sum_log_ratios(ll_after_block, ilr2) 
      }
    }
    psis_part <- suppressWarnings(psis(logratio))
    k <- pareto_k_values(psis_part)
    ks[i] <- k
    if (k > k_thres) {
      # refit the model based on the first i observations
      i_star <- i
      refits <- c(refits, i)
      ind_rm <- (i + 1):min(i + B, N)
      newdf <- df
      if (max(ind_rm) < N) {
        # responses of the left-out block need to be coded as 
        # missing to retain a single consistent time series
        newdf$y[ind_rm] <- NA
      } else {
        newdf <- newdf[-ind_rm, , drop = FALSE]
      }
      fit_star <- update(fit, newdata = newdf, recompile = FALSE, refresh = 0)
      conv[[i]] <- convergence_summary(fit_star)
      # perform exact LFO for the ith observation
      ll <- log_lik(fit_star, newdata = df[ioos, , drop = FALSE], oos = oos)
      loglik[, i + 1] <- ll[, i + 1]
      out[i] <- lfo_criterion(
        fit_star, data = df[ioos, , drop = FALSE], 
        oos = oos, criterion = criterion
      )
    } else {
      # PSIS approximate LFO is possible
      out[i] <- lfo_criterion(
        fit_star, data = df[ioos, , drop = FALSE], oos = oos,
        criterion = criterion, psis = psis_part
      )
    }
  }
  attr(out, "ks") <- ks
  attr(out, "refits") <- refits
  attr(out, "conv") <- conv
  out
}

compute_lfo <- function(fit, type = c("exact", "approx"), file = NULL, ...) {
  type <- match.arg(type)
  if (!is.null(file) && file.exists(file)) {
    return(readRDS(file))
  }
  lfo_fun <- get(paste0(type, "_lfo"), mode = "function")
  out <- lfo_fun(fit, ...)
  if (!is.null(file)) {
    saveRDS(out, file)
  }
  out
}

lfo_criterion <- function(fit, data, oos, criterion = c("elpd", "rmse"),
                          psis = NULL, ...) {
  stopifnot(is.brmsfit(fit))
  data <- as.data.frame(data)
  stopifnot(all(oos %in% seq_len(NROW(data))))
  criterion <- match.arg(criterion)
  stopifnot(is.null(psis) || loo:::is.psis(psis))
  if (criterion == "elpd") {
    ll <- log_lik(fit, newdata = data, oos = oos, ...)
    sum_ll <- rowSums(ll[, oos, drop = FALSE])
    if (is.null(psis)) {
      out <- log_mean_exp(sum_ll)
    } else {
      lw <- weights(psis, normalize = TRUE)[, 1]
      out <- log_sum_exp(lw + sum_ll)
    }
  } else if (criterion == "rmse") {
    pp <- posterior_predict(fit, newdata = data, oos = oos, ...)
    y <- brms:::get_y(fit, newdata = data)
    if (is.null(psis)) {
      out <- rmse(pp[, oos, drop = FALSE], y[oos])
    } else {
      w <- weights(psis, log = FALSE, normalize = TRUE)[, 1]
      out <- rmse(pp[, oos, drop = FALSE], y[oos], weights = w)
    }
  }
  out
}

convergence_summary <- function(x) {
  require(dplyr)
  list(
    nuts_params = nuts_params(x) %>%
      spread("Parameter", "Value") %>%
      group_by(Chain) %>%
      summarise(
        n_divergent = sum(divergent__),
        max_treedepth = max(treedepth__),
        stepsize = mean(stepsize__),
        mean_n_leapfrog = mean(n_leapfrog__)
      ),
    max_rhat = max(rhat(x)),
    min_neff_ratio = min(neff_ratio(x))
  )
}

fit_model <- function(model, N, ...) {
  require(brms)
  time <- seq_len(N)
  stime <- scale_unit_interval(time)
  df <- data.frame(time = time, stime = stime)
  ar_prior <- prior(normal(0, 0.5), class = "ar")
  ar_autocor <- cor_ar(~ stime, p = 2)
  if (model == "constant") {
    df$y <- rnorm(N)
    fit <- brm(y | mi() ~ 1, data = df, refresh = 0, ...)
  } else if (model == "linear") {
    df$y <- 17 * stime + rnorm(N) 
    fit <- brm(y | mi() ~ stime, data = df, refresh = 0, ...)
  } else if (model == "quadratic") {
    df$y <- 17 * stime - 25 * stime^2 + rnorm(N)
    fit <- brm(y | mi() ~ stime + I(stime^2), data = df, refresh = 0, ...)
  } else if (model == "AR2-only") {
    df$y <- as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y | mi() ~ 1, data = df, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else if (model == "AR2-linear") {
    df$y <- 17 * stime +
      as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y | mi() ~ stime, data = df, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else if (model == "AR2-quadratic") {
    df$y <- 17 * stime - 25 * stime^2 + 
      as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y | mi() ~ stime + I(stime^2), 
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
  N <- conditions$N[j]
  M <- conditions$M[j]
  L <- conditions$L[j]
  B <- conditions$B[j]
  k_thres <- conditions$k_thres[j]
  model <- conditions$model[j]
  fit <- fit_model(model = model, N = N, ...)

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
