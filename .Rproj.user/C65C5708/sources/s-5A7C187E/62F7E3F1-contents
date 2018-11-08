source("sim_functions.R")
library(tidyverse)
library(foreach)
library(doParallel)
chains <- 1
iter <- 5000
warmup <- 1000
control <- list(adapt_delta = 0.95)
# SEED <- 121
# set.seed(SEED)

nsims <- 1
conditions <- expand.grid(
  N = 200, M = 1, L = 25, B = NA,
  model = c("constant", "AR2_only", "AR2_linear", "AR2_quadratic"),
  k_thres = c(0.5, 0.6, 0.7),
  sim = seq_len(nsims)
)
conditions$res <- list(list())

J <- seq_len(nrow(conditions))

cl <- makeCluster(1)
registerDoParallel(cl)

conditions$res[I] <- 
  foreach(j = J, .packages = c("brms", "loo")) %dopar% 
  sim_fun(j, conditions, chains = chains, iter = iter, 
          warmup = warmup, control = control)

stopCluster(cl)

file <- "results/lfo_models.rds"
if (file.exists(file)) {
  old_conditions <- read_rds(file)
  conditions <- bind_rows(old_conditions, conditions)
}
write_rds(conditions, file)
