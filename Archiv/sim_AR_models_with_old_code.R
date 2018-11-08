source("sim_helpers.R")
library(tidyverse)
library(loo)
library(brms)
CHAINS <- 1
ITER <- 2000
control <- list(adapt_delta = 0.95)
SEED <- 12
set.seed(SEED)

N <- 200
time <- 1:200
stime <- scale_unit_interval(time)
L <- 25
k_thres <- 0.6

nsims <- 20
res <- vector("list", nsims)
for (j in seq_along(res)) {
  # simulate from an AR(2) model
  y <- arima.sim(list(ar = c(0.5, 0.3)), N) - 25 * stime^2 + 17 * stime
  df <- data.frame(y = as.numeric(y), time = time, stime = stime)
  
  # ggplot(df, aes(x = time, y = y)) + 
  #   geom_point(size = 1)
  
  fit <- brm(
    y ~ 1 + stime + (stime^2), 
    data = df, 
    autocor = cor_ar(~stime, p = 2), 
    prior = prior(normal(0, 0.5), class = "ar"),
    chains = CHAINS, iter = ITER,
    control = control
  )
  
  # preds <- posterior_predict(fit)
  # preds <- cbind(
  #   Estimate = colMeans(preds), 
  #   Q5 = apply(preds, 2, quantile, probs = 0.05),
  #   Q95 = apply(preds, 2, quantile, probs = 0.95)
  # )
  # ggplot(cbind(df, preds), aes(x = time, y = Estimate)) +
  #   geom_smooth(aes(ymin = Q5, ymax = Q95), stat = "identity", size = 0.5) +
  #   geom_point(aes(y = y))
  
  # LOO-CV
  loo_cv <- loo(log_lik(fit)[, (L + 1):N])
  
  # Exact LFO-CV
  exact_elpds_1sap <- exact_lfo(fit, M = 1, L = L, k_thres = k_thres)
  exact_elpd_1sap <- summarize_elpds(exact_elpds_1sap)
  
  # loglik_exact <- matrix(nrow = nsamples(fit), ncol = N)
  # for (i in N:max(L + 1, 2)) {
  #   fit_i <- update(fit, newdata = df[-(i:N), ], recompile = FALSE)
  #   loglik_exact[, i] <- log_lik(fit_i, newdata = df[1:i, ])[, i]
  # }
  # exact_elpds_1sap <- na.omit(apply(loglik_exact, 2, log_mean_exp))
  # exact_elpd_1sap <- c(
  #   ELPD = sum(exact_elpds_1sap),
  #   SE = sqrt(length(exact_elpds_1sap) * var(exact_elpds_1sap))
  # )
  
  # Approximate LFO-CV
  approx_elpds_1sap <- approx_lfo(fit, M = 1, L = L, k_thres = k_thres)
  approx_elpd_1sap <- summarize_elpds(approx_elpds_1sap)
  
  # loglik <- matrix(nrow = nsamples(fit), ncol = N)
  # approx_elpds_1sap <- rep(NA, N)
  # fit_part <- fit
  # ids <- N:(L + 1)
  # i_refit <- N
  # refits <- NULL
  # ks <- NULL
  # for (i in ids) {
  #   loglik[, i] <- log_lik(fit_part)[, i]
  #   logratio <- sum_log_ratios(loglik, i:i_refit)
  #   psis_part <- suppressWarnings(psis(logratio))
  #   k <- pareto_k_values(psis_part)
  #   ks <- c(ks, k)
  #   if (k > k_thres) {
  #     # refit the model based on the first i-1 observations
  #     i_refit <- i
  #     refits <- c(refits, i)
  #     fit_part <- update(fit_part, newdata = df[1:(i - 1), ], recompile = FALSE)
  #     loglik[, i] <- log_lik(fit_part, newdata = df[1:i, ])[, i]
  #     approx_elpds_1sap[i] <- log_mean_exp(loglik[, i])
  #   } else {
  #     lw_i <- weights(psis_part, normalize = TRUE)[, 1]
  #     approx_elpds_1sap[i] <- log_sum_exp(lw_i + loglik[, i])
  #   }
  # } 
  # approx_elpds_1sap <- na.omit(approx_elpds_1sap)
  # approx_elpd_1sap <- c(
  #   ELPD = sum(approx_elpds_1sap),
  #   SE = sqrt(length(approx_elpds_1sap) * var(approx_elpds_1sap))
  # )
  
  # rbind_print(
  #   "LOO" = loo_cv$estimates["elpd_loo", ],
  #   "LFO" = exact_elpd_1sap
  # )
  
  # dat_elpd <- data.frame(
  #   loo_elpd = loo_cv$pointwise[, 1],
  #   exact_elpd = exact_elpds_1sap,
  #   time = seq_along(exact_elpds_1sap)
  # )
  # 
  # ggplot(dat_elpd, aes(x = loo_elpd, y = exact_elpd)) +
  #   geom_abline(color = "gray30") +
  #   geom_point(size = 2) +
  #   labs(x = "LOO ELPDs", y = "Exact ELPDs")
  # 
  # dat_elpd$diff <- dat_elpd$exact_elpd - dat_elpd$loo_elpd
  # ggplot(dat_elpd, aes(x = time, y = diff)) +
  #   geom_point(size = 2) +
  #   labs(x = "time", y = "ELPD Diffs")
  
  res[[j]] <- list(
    df = df, fit = fit, loo_cv = loo_cv,
    exact_elpds_1sap = exact_elpds_1sap, 
    exact_elpd_1sap = exact_elpd_1sap,
    approx_elpds_1sap = approx_elpds_1sap, 
    approx_elpd_1sap = approx_elpd_1sap
  )
}

write_rds(res, "results/ar_models_1sap.rds")
