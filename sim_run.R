source("sim_functions.R")
library(dplyr)
library(tidyr)
library(readr)
library(foreach)
library(doParallel)
chains <- 1
iter <- 5000
warmup <- 1000
control <- list(adapt_delta = 0.95)

nsims <- 100
conditions <- expand.grid(
  N = 200, 
  M = c(1, 4),
  L = 25, 
  B = NA,
  model = c(
    "constant", "linear", "quadratic",
    "AR2-only", "AR2-linear", "AR2-quadratic"
  ),
  k_thres = c(0.5, 0.6, 0.7),
  sim = seq_len(nsims)
)
conditions$res <- list(list())

cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)

J <- seq_len(nrow(conditions))
conditions$res[J] <- 
  foreach(j = J, .packages = c("brms", "loo")) %dopar% 
  sim_fun(j, conditions, chains = chains, iter = iter, 
          warmup = warmup, control = control)

stopCluster(cl)

file <- "results/lfo_sims.rds"
if (file.exists(file)) {
  old_conditions <- read_rds(file)
  conditions <- bind_rows(old_conditions, conditions)
}
write_rds(conditions, file)
