# some helper functions we'll use throughout
SW <- suppressWarnings
SM <- suppressMessages

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
sum_log_ratios <- function(ll, ids = NULL, mode) {
  if (!is.null(ids)) {
    ll <- ll[, ids, drop = FALSE]
  }
  out <- rowSums(ll)
  if (mode == "backward") {
    out <- -out
  }
  out
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

plot_ks <- function(ks, ids, k_thres = 0.7) {
  dat_ks <- data.frame(ks = ks, ids = ids)
  ggplot(dat_ks, aes(x = ids, y = ks)) + 
    geom_point(aes(color = ks > k_thres), shape = 3, show.legend = FALSE) + 
    geom_hline(yintercept = k_thres, linetype = 2, color = "red2") + 
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
# @param fit a brmsfit object
# @param M steps to predict into the future
# @param L minimal number of observations required for fitting the model
# @return A vector of pointwise ELPD values
exact_lfo <- function(fit, M, L, criterion = c("elpd", "rmse"), 
                      factorize = FALSE, ...) {
  require(brms)
  stopifnot(is.brmsfit(fit))
  criterion <- match.arg(criterion)
  data <- model.frame(fit)
  N <- NROW(data)
  
  # compute exact LFO likelihoods
  out <- rep(NA, N)
  conv <- vector("list", N)
  ids <- L:(N - M)
  for (i in ids) {
    tmp <- lfo_step_pred(
      i = i, fit_star = fit, data = data, M = M, 
      criterion = criterion, factorize = factorize
    )
    out[i] <- tmp$crit
  }
  attr(out, "conv") <- conv
  out
}

# Perform approximate leave-future-out cross-validation
# @param fit a brmsfit object
# @param M steps to predict into the future
# @param L minimal number of observations required for fitting the model
# @return A vector of pointwise ELPD values
approx_lfo <- function(fit, M, L, k_thres = 0.7, 
                       mode = c("forward", "backward", "combined"),
                       criterion = c("elpd", "rmse"),
                       factorize = FALSE, ...) {
  require(brms)
  require(loo)
  stopifnot(is.brmsfit(fit))
  mode <- match.arg(mode)
  criterion <- match.arg(criterion)
  data <- model.frame(fit)
  N <- NROW(data)
  
  # compute approximate LFO likelihoods
  out <- ks <- reffs <- rep(NA, N)
  conv <- vector("list", N)
  refits <- numeric(0)
  if (mode %in% c("forward", "combined")) {
    # move forward in time
    # need to refit the model with the first L observations
    reffs[L] <- 1
    tmp <- lfo_step_pred(
      i = L, fit_star = fit, data = data, M = M, 
      criterion = criterion, factorize = factorize
    )
    fit_star <- tmp$fit_star
    i_star <- tmp$i_star
    ks[L] <- tmp$k
    refits[L] <- tmp$refit
    conv[[L]] <- tmp$conv
    out[L] <- tmp$crit 
    # start from L + 1 as we already handled L above
    ids <- (L + 1):(N - M)
  } else {
    # move backward in time
    fit_star <- fit
    i_star <- N
    ids <- (N - M):L
  }

  for (i in ids) {
    if (mode == "combined") {
      i_star_old <- i_star
    }
    psis_i <- lfo_step_psis(
      i = i, fit_star = fit_star, i_star = i_star, 
      data = data, mode = mode, factorize = factorize
    )
    reffs[i] <- attr(psis_i, "r_eff_lr")
    tmp <- lfo_step_pred(
      i = i, fit_star = fit_star, i_star = i_star, 
      data = data, M = M, psis = psis_i, k_thres = k_thres,
      criterion = criterion, factorize = factorize
    )
    fit_star <- tmp$fit_star
    i_star <- tmp$i_star
    ks[i] <- tmp$k
    refits[i] <- tmp$refit
    conv[[i]] <- tmp$conv
    out[i] <- tmp$crit
    
    if (i_star == i && mode == "combined") {
      # model was just refitted now also move backwards in time
      ids_bw <- (i_star - 1):(i_star_old + 1) 
      for (j in ids_bw) {
        psis_j <- lfo_step_psis(
          i = j, fit_star = fit_star, i_star = i_star, 
          data = data, mode = "backward", factorize = factorize
        )
        k_j <- pareto_k_values(psis_j)
        if (k_j > k_thres) {
          # stop backward in forward if a refit is necessary
          break
        } else {
          tmp <- lfo_step_pred(
            i = j, fit_star = fit_star, i_star = i_star, 
            data = data, M = M, psis = psis_j, k_thres = k_thres,
            criterion = criterion, factorize = factorize
          )
          # TODO: use multiple proposal importance sampling
          out[j] <- (out[j] + tmp$crit) / 2 
        }
      }
    }
  }
  attr(out, "ks") <- ks
  attr(out, "reffs") <- reffs
  attr(out, "refits") <- refits
  attr(out, "conv") <- conv
  out
}

lfo_step_psis <- function(i, fit_star, data, i_star, mode,
                          factorize = FALSE) {
  # compute psis over observations in J_i 
  if (mode %in% c("forward", "combined")) {
    J_i <- (i_star + 1):i 
  } else {
    J_i <- (i + 1):i_star
  }
  if (factorize) {
    data_psis <- data[J_i, , drop = FALSE]
    J_i <- seq_len(NROW(data_psis))
  } else {
    data_psis <- data[1:max(J_i), , drop = FALSE]
  }

  ll_psis <- log_lik(fit_star, newdata = data_psis)
  logratio <- sum_log_ratios(ll_psis, ids = J_i, mode = mode)
  
  # compute relative efficiency of logratios
  chains <- fit_star$fit@sim$chains
  chain_id <- rep(seq_len(chains), each = NROW(logratio) / chains)
  r_eff_lr <- loo::relative_eff(logratio, chain_id = chain_id)
  
  r_eff <- loo::relative_eff(1 / exp(logratio), chain_id = chain_id)
  psis <- SW(psis(logratio, r_eff = r_eff))
  attr(psis, "r_eff_lr") <- r_eff_lr
  psis
}

lfo_step_pred <- function(i, fit_star, data, M, criterion, factorize = FALSE,
                          k_thres = Inf, i_star = i, psis = NULL, ...) {
  print(i)
  # observations to predict
  oos <- (i + 1):(i + M)
  data_pred <- data[1:(i + M), , drop = FALSE]
  
  if (!is.null(psis)) {
    k <- pareto_k_values(psis) 
  } else {
    k <- 0
  }
  if (is.null(psis) || k > k_thres) {
    # refit the model based on the first i observations
    refit <- TRUE
    i_star <- i
    fit_star <- SM(update(
      fit_star, newdata = data[1:i, , drop = FALSE], 
      recompile = FALSE, refresh = 0, ...
    ))
    conv <- convergence_summary(fit_star)
    # perform exact LFO-CV for the ith observation
    crit <- comp_criterion(
      fit_star, data = data_pred, oos = oos, 
      criterion = criterion, factorize = factorize
    )
  } else {
    # PSIS approximate LFO-CV is possible
    refit <- FALSE
    conv <- list()
    crit <- comp_criterion(
      fit_star, data = data_pred, oos = oos,
      criterion = criterion, psis = psis,
      factorize = factorize
    )
  }
  list(
    crit = crit,
    i_star = i_star,
    fit_star = fit_star,
    k = k, 
    conv = conv,
    refit = refit
  )
}

comp_criterion <- function(fit, data, oos, criterion = c("elpd", "rmse"),
                           factorize = FALSE, psis = NULL, ...) {
  stopifnot(is.brmsfit(fit))
  data <- as.data.frame(data)
  stopifnot(all(oos %in% seq_len(NROW(data))))
  criterion <- match.arg(criterion)
  stopifnot(is.null(psis) || loo:::is.psis(psis))
  if (factorize) {
    # we don't need past observations to make predictions
    data <- data[oos, , drop = FALSE]
    oos <- seq_len(NROW(data))
  }
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
  require(tidyr)
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

fit_model <- function(model, N, ...) {
  require(brms)
  time <- seq_len(N)
  stime <- scale_unit_interval(time)
  data <- data.frame(time = time, stime = stime)
  ar_prior <- prior(normal(0, 0.5), class = "ar")
  ar_autocor <- cor_ar(~ stime, p = 2)
  if (model == "constant") {
    data$y <- rnorm(N)
    fit <- brm(y | mi() ~ 1, data = data, refresh = 0, ...)
  } else if (model == "linear") {
    data$y <- 17 * stime + rnorm(N) 
    fit <- brm(y | mi() ~ stime, data = data, refresh = 0, ...)
  } else if (model == "quadratic") {
    data$y <- 17 * stime - 25 * stime^2 + rnorm(N)
    fit <- brm(y | mi() ~ stime + I(stime^2), data = data, refresh = 0, ...)
  } else if (model == "AR2-only") {
    data$y <- as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y | mi() ~ 1, data = data, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else if (model == "AR2-linear") {
    data$y <- 17 * stime +
      as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y | mi() ~ stime, data = data, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else if (model == "AR2-quadratic") {
    data$y <- 17 * stime - 25 * stime^2 + 
      as.numeric(arima.sim(list(ar = c(0.5, 0.3)), N))
    fit <- brm(
      y | mi() ~ stime + I(stime^2), 
      data = data, prior = ar_prior,
      autocor = ar_autocor, refresh = 0, ...
    )
  } else {
    stop("Model '", model, "' is not supported.")
  }
  fit
}

sim_fun <- function(j, conditions, ...) {
  # simulate data and fit the corresponding model
  message("Simulating LFO-CV for condition ", j)
  N <- conditions$N[j]
  M <- conditions$M[j]
  L <- conditions$L[j]
  k_thres <- conditions$k_thres[j]
  model <- conditions$model[j]
  fit <- fit_model(model = model, N = N, ...)
  loo_cv <- loo(log_lik(fit)[, (L + 1):N])
  
  lfo_exact_elpds <- exact_lfo(fit, M = M, L = L)
  lfo_approx_fw_elpds <- approx_lfo(
    fit, M = M, L = L, k_thres = k_thres, mode = "forward"
  )
  lfo_approx_bw_elpds <- approx_lfo(
    fit, M = M, L = L, k_thres = k_thres, mode = "backward"
  )
  # the combined mode is not yet fully implemented
  # lfo_approx_cb_elpds <- approx_lfo(
  #   fit, M = M, L = L, k_thres = k_thres, mode = "combined"
  # )
  
  # return all relevant information
  list(
    # storing 'fit' requires too much space
    loo_cv = loo_cv,
    lfo_exact_elpds = lfo_exact_elpds,
    lfo_approx_fw_elpds = lfo_approx_fw_elpds,
    lfo_approx_bw_elpds = lfo_approx_bw_elpds
    # lfo_approx_cb_elpds = lfo_approx_cb_elpds
  )
}

sim_fun_rmse <- function(j, conditions, ...) {
  # simulate data and fit the corresponding model
  message("Simulating LFO-CV for condition ", j)
  N <- conditions$N[j]
  M <- conditions$M[j]
  L <- conditions$L[j]
  k_thres <- conditions$k_thres[j]
  model <- conditions$model[j]
  fit <- fit_model(model = model, N = N, ...)
  
  lfo_exact_rmses <- exact_lfo(fit, M = M, L = L, criterion = "rmse")
  lfo_approx_fw_rmses <- approx_lfo(
    fit, M = M, L = L, k_thres = k_thres,
    mode = "forward", criterion = "rmse"
  )
  
  # return all relevant information
  list(
    # storing 'fit' requires too much space
    lfo_exact_rmses = lfo_exact_rmses,
    lfo_approx_fw_rmses = lfo_approx_fw_rmses
  )
}
