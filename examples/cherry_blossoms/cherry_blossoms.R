source("sim_functions.R")

cherry <- read.csv("examples/cherry_blossoms/cherry_blossoms.csv")
cherry_temp <- cherry[!is.na(cherry$temp), ]
cherry_doy <- cherry[!is.na(cherry$doy), ]

library(brms)
options(mc.cores = 2)

fit_cb <- brm(
  formula = bf(doy ~ s(year, k = 40)),
  data = cherry_doy, chain = 2, seed = 583829, 
  control = list(adapt_delta = 0.99)
)
plot(marginal_effects(fit_cb), points = TRUE)

L <- 100
M <- 1

# perform exact LFO
exact_lfo_cb <- exact_lfo(fit_cb, M = M, L = L)
saveRDS(exact_lfo_cb, file = "examples/exact_lfo_cb.rds")
summarize_elpds(exact_lfo_cb)

# perform approximate LFO
approx_lfo_cb <- approx_lfo(fit_cb, M = M, L = L)
saveRDS(approx_lfo_cb, file = "examples/approx_lfo_cb.rds")
summarize_elpds(approx_lfo_cb)

# perform approximate LOO
(loo_cb <- loo(fit_cb, newdata = cherry_doy[-seq_len(L), ]))

