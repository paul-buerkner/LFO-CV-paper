library(tidyverse)
source("sim_functions.R")
colors <- brms:::viridis6()[c(6, 2, 3, 4)]
theme_set(theme_bw())
# colors <- unname(unlist(bayesplot::color_scheme_get()[c(6, 4, 2)]))
# theme_set(bayesplot::theme_default(base_family = "sans"))

# helper function to extract the number of refits
nrefits <- function(x) {
  sum(attr(x, "refits"), na.rm = TRUE)
}

# data preparation
mlevels <- c(
  "constant", "linear", "quadratic",
  "AR2-only", "AR2-linear", "AR2-quadratic"
)

lfo_sims <- read_rds("results/lfo_sims.rds") %>%
  as_tibble() %>%
  filter(lengths(res) > 0) %>%
  mutate(
    model = factor(model, levels = mlevels),
    k_thres = paste0("k = ", k_thres),
    # compute sum of pointwise ELPD values
    elpd_loo = map_dbl(res, ~ .$loo_cv$estimates["elpd_loo", 1]),
    elpd_exact_lfo = map_dbl(res, ~ summarize_elpds(.$lfo_exact_elpds)[1]),
    elpd_approx_fw_lfo = map_dbl(res, ~ summarize_elpds(.$lfo_approx_fw_elpds)[1]),
    elpd_approx_bw_lfo = map_dbl(res, ~ summarize_elpds(.$lfo_approx_bw_elpds)[1]),
    elpd_approx_cb_lfo  = map_dbl(res, ~ summarize_elpds(.$lfo_approx_cb_elpds)[1]),
    # compute differences to the exact values
    elpd_diff_bw_lfo = elpd_approx_bw_lfo - elpd_exact_lfo,
    elpd_diff_fw_lfo = elpd_approx_fw_lfo - elpd_exact_lfo,
    elpd_diff_cb_lfo = elpd_approx_cb_lfo - elpd_exact_lfo,
    elpd_diff_loo = elpd_loo - elpd_exact_lfo,
    nrefits_bw = map_dbl(res, ~ nrefits(.$lfo_approx_bw_elpds)),
    nrefits_fw = map_dbl(res, ~ nrefits(.$lfo_approx_fw_elpds)),
    nrefits_cb = map_dbl(res, ~ nrefits(.$lfo_approx_cb_elpds)),
    rel_nrefits_bw = nrefits_bw / (N - L),
    rel_nrefits_fw = nrefits_fw / (N - L),
    rel_nrefits_cb = nrefits_cb / (N - L)
  ) %>%
  select(-res)


# LFO-1SAP plots ----------------------------------------------------------
lfo_res_1sap <- lfo_sims %>% filter(is.na(B), M == 1)

lfo_res_1sap %>%
  group_by(model, k_thres) %>%
  summarise(
    mean_diff_bw_lfo = round(mean(elpd_diff_bw_lfo), 2), 
    mean_diff_fw_lfo = round(mean(elpd_diff_fw_lfo), 2), 
    mean_diff_cb_lfo = round(mean(elpd_diff_cb_lfo), 2), 
    mean_diff_loo = round(mean(elpd_diff_loo), 2),
    sd_diff_bw_lfo = round(sd(elpd_diff_bw_lfo), 2),
    sd_diff_fw_lfo = round(sd(elpd_diff_fw_lfo), 2),
    sd_diff_cb_lfo = round(sd(elpd_diff_cb_lfo), 2),
    sd_diff_loo = round(sd(elpd_diff_loo), 2)
  )

# plot differences to exact LFO
lfo_res_1sap %>% 
  select(starts_with("elpd_diff"), model, k_thres, M) %>%
  gather("Type", "elpd_diff", starts_with("elpd_diff")) %>%
  ggplot(aes(x = elpd_diff, y = ..density.., fill = Type)) +
  geom_histogram(alpha = 0.7) +
  facet_grid(model ~ k_thres, scales = "free_y") +
  scale_fill_manual(
    values = colors,
    labels = c(
      "PSIS-LFO-CV (backward)",
      "PSIS-LFO-CV (combined)",
      "PSIS-LFO-CV (forward)",
      "PSIS-LOO-CV"
    )
  ) +
  labs(x = 'Difference to exact LFO-CV', y = "Density") +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(legend.position = "bottom") +
  NULL

ggsave("plots/LFO_1SAP_AR_models_ELPD.jpeg", width = 7, height = 6)


# plot relative number of refits
lfo_res_1sap %>% 
  select(starts_with("rel_nrefits"), model, k_thres) %>%
  gather("Type", "rel_nrefits", starts_with("rel_nrefits")) %>%
  ggplot(aes(x = rel_nrefits, y = ..density.., fill = Type)) +
  facet_grid(model ~ k_thres) + 
  geom_vline(
    aes(xintercept = rel_nrefits),
    linetype = 2,
    data = lfo_res %>%
      select(starts_with("rel_nrefits"), model, k_thres) %>%
      gather("Type", "rel_nrefits", starts_with("rel_nrefits")) %>%
      group_by(model, k_thres, Type) %>%
      summarise(rel_nrefits = mean(rel_nrefits))
  ) +
  geom_histogram(alpha = 0.7) +
  scale_fill_manual(
    values = colors,
    labels = c(
      "PSIS-LFO-CV (backward)",
      "PSIS-LFO-CV (combined)",
      "PSIS-LFO-CV (forward)",
      "PSIS-LOO-CV"
    )
  ) +
  xlab("Relative number of refits based on N = 200 and L = 25") +
  ylab("Count") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("plots/LFO_1SAP_AR_models_rel_refits.jpeg", width = 6, height = 6)


# LFO-4SAP plots ---------------------------------------------------------
lfo_res_4sap <- lfo_sims %>% 
  filter(is.na(B), M == 4)

# plot differences to exact LFO
lfo_res_4sap %>% 
  select(starts_with("elpd_diff"), model, k_thres, M) %>%
  gather("Type", "elpd_diff", starts_with("elpd_diff")) %>%
  filter(Type != "elpd_diff_loo") %>%
  ggplot(aes(x = elpd_diff, y = ..density.., fill = Type)) +
  geom_histogram(alpha = 0.7) +
  facet_grid(model ~ k_thres, scales = "free_y") +
  scale_fill_manual(
    values = colors,
    labels = c(
      "PSIS-LFO-CV (backward)",
      "PSIS-LFO-CV (combined)",
      "PSIS-LFO-CV (forward)",
      "PSIS-LOO-CV"
    )
  ) +
  labs(x = 'Difference to exact LFO-CV', y = "Density") +
  geom_vline(xintercept = 0, linetype = 2) +
  xlim(c(-20, 50)) +
  theme(legend.position = "bottom") +
  NULL

ggsave("plots/LFO_4SAP_AR_models_ELPD.jpeg", width = 7, height = 6)


# plot relative number of refits
lfo_res_4sap %>% 
  select(starts_with("rel_nrefits"), model, k_thres) %>%
  gather("Type", "rel_nrefits", starts_with("rel_nrefits")) %>%
  ggplot(aes(x = rel_nrefits, y = ..density.., fill = Type)) +
  facet_grid(model ~ k_thres) + 
  geom_vline(
    aes(xintercept = rel_nrefits),
    linetype = 2,
    data = lfo_res %>%
      select(starts_with("rel_nrefits"), model, k_thres) %>%
      gather("Type", "rel_nrefits", starts_with("rel_nrefits")) %>%
      group_by(model, k_thres, Type) %>%
      summarise(rel_nrefits = mean(rel_nrefits))
  ) +
  geom_histogram(alpha = 0.7) +
  scale_fill_manual(
    values = colors,
    labels = c(
      "PSIS-LFO-CV (backward)",
      "PSIS-LFO-CV (combined)",
      "PSIS-LFO-CV (forward)",
      "PSIS-LOO-CV"
    )
  ) +
  xlab("Relative number of refits based on N = 200 and L = 25") +
  ylab("Count") +
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("plots/LFO_4SAP_AR_models_rel_refits.jpeg", width = 6, height = 6)

