#library
library(tidybayes)
library(ggplot2)
library(patchwork)
library(dplyr)


# ── Biting models ──────────────────────────────────────────────────

bite_dist_brms <- brm(
  Seafloor_bites_binary ~ avg_distance_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = full_data_strip,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),    class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),  class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "bite_dist_brms"  # saves model so you don't refit
)

summary(bite_dist_brms)

bite_height_brms <- brm(
  Seafloor_bites_binary ~ avg_norm_height_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = full_data_strip,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "bite_height_brms"
)

summary(bite_height_brms)

bite_hitrate_brms <- brm(
  Seafloor_bites_binary ~ hit_rate_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = full_data_strip,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "bite_hitrate_brms"
)
summary(bite_hitrate_brms)

# ── Foraging initiation models ─────────────────────────────────────

f_init_dist_brms <- brm(
  next_state_different ~ next_avg_distance_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_foraging_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "finit_dist_brms"
)
summary(f_init_dist_brms)

f_init_height_brms <- brm(
  next_state_different ~ next_avg_norm_height_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_foraging_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "finit_height_brms"
)
summary(f_init_height_brms)

f_init_hitrate_brms <- brm(
  next_state_different ~ next_hit_rate_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_foraging_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "finit_hitrate_brms"
)

summary(f_init_hitrate_brms)

# ── Traversing initiation models ─────────────────────────────────────
t_init_dist_brms <- brm(
  next_state_different ~ next_avg_distance_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_traverse_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "tinit_dist_brms"
)
summary(t_init_dist_brms)

t_init_height_brms <- brm(
  next_state_different ~ next_avg_norm_height_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_traverse_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "tinit_height_brms"
)
summary(t_init_height_brms)

t_init_hitrate_brms <- brm(
  next_state_different ~ next_hit_rate_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_traverse_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "tinit_hitrate_brms"
)
summary(t_init_hitrate_brms )

###############
#plot

# Recover scale attributes for biting data
attr_dist   <- list(mean = mean(full_data_strip$avg_distance,    na.rm = TRUE),
                    sd   = sd(full_data_strip$avg_distance,      na.rm = TRUE))
attr_height <- list(mean = mean(full_data_strip$avg_norm_height, na.rm = TRUE),
                    sd   = sd(full_data_strip$avg_norm_height,   na.rm = TRUE))
attr_hit    <- list(mean = mean(full_data_strip$hit_rate,        na.rm = TRUE),
                    sd   = sd(full_data_strip$hit_rate,          na.rm = TRUE))

# ── Generic plot function ──────────────────────────────────────────

plot_brms_fit <- function(model, data,
                          x_var_z,
                          x_var_raw,
                          x_mean, x_sd,
                          y_var,
                          x_label, y_label,
                          social_var,
                          ndraws = 500,
                          n_bins = 20) {
  
  # ── Raw counts at 0 and 1 ────────────────────────────────────────
  raw_counts <- data %>%
    mutate(x_bin = cut(.data[[x_var_raw]], breaks = n_bins)) %>%
    group_by(x_bin, subject, .data[[social_var]], .data[[y_var]]) %>%
    summarise(
      x_mid = mean(.data[[x_var_raw]], na.rm = TRUE),
      n     = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(x_mid)) %>%
    # Offset social and asocial vertically so both visible
    mutate(
      y_offset = case_when(
        .data[[social_var]] == "Asocial" & .data[[y_var]] == 0 ~ -0.02,
        .data[[social_var]] == "Asocial" & .data[[y_var]] == 1 ~  0.98,
        .data[[social_var]] == "Social"  & .data[[y_var]] == 0 ~  0.02,
        .data[[social_var]] == "Social"  & .data[[y_var]] == 1 ~  1.02
      )
    )
  
  # ── Prediction grid ──────────────────────────────────────────────
  # ── Per-combination data ranges ───────────────────────────────────
  x_ranges <- data %>%
    group_by(subject, .data[[social_var]]) %>%
    summarise(
      x_min = min(.data[[x_var_raw]], na.rm = TRUE),
      x_max = max(.data[[x_var_raw]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # ── Prediction grid — bounded per combination ─────────────────────
  pred_grid <- x_ranges %>%
    group_by(subject, .data[[social_var]]) %>%
    tidyr::expand(x_raw = seq(x_min[1], x_max[1], length.out = 100)) %>%
    ungroup() %>%
    mutate(!!x_var_z := (x_raw - x_mean) / x_sd)
  
  if ("time_in_state" %in% colnames(data)) {
    pred_grid$time_in_state <- median(data$time_in_state, na.rm = TRUE)
  }
  
  names(pred_grid)[names(pred_grid) == "social"] <- social_var
  
  pred_grid$Trial_ID    <- NA
  pred_grid$site        <- NA
  pred_grid$folder_date <- NA
  
  # ── Posterior predictions ────────────────────────────────────────
  pred_draws <- pred_grid %>%
    add_epred_draws(model, re_formula = NA, ndraws = ndraws)
  
  pred_summary <- pred_draws %>%
    group_by(x_raw, subject, .data[[social_var]]) %>%
    summarise(
      median = median(.epred),
      lower  = quantile(.epred, 0.025),
      upper  = quantile(.epred, 0.975),
      .groups = "drop"
    )
  
  # ── Plot ─────────────────────────────────────────────────────────
  ggplot() +
    geom_point(
      data = raw_counts,
      aes(x     = x_mid,
          y     = y_offset,
          color = .data[[social_var]],
          size  = n),
      alpha = 0.4
    ) +
    geom_ribbon(
      data = pred_summary,
      aes(x    = x_raw,
          ymin = lower,
          ymax = upper,
          fill = .data[[social_var]]),
      alpha = 0.2
    ) +
    geom_line(
      data = pred_summary,
      aes(x     = x_raw,
          y     = median,
          color = .data[[social_var]]),
      linewidth = 1
    ) +
    scale_size_continuous(
      name  = "N",
      range = c(0.5, 6)
    ) +
    facet_wrap(~ subject) +
    scale_color_manual(
      values = c("Asocial" = "#F8766D",
                 "Social"  = "#00BFC4")
    ) +
    scale_fill_manual(
      values = c("Asocial" = "#F8766D",
                 "Social"  = "#00BFC4")
    ) +
    labs(x     = x_label,
         y     = y_label,
         color = "Social state",
         fill  = "Social state",
         size  = "N") +
    theme_classic(base_size = 18) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          legend.position  = "none")
}

# ── Scale attributes ───────────────────────────────────────────────

# Biting — recovered manually
b_dist_mean   <- mean(full_data_strip$avg_distance,    na.rm = TRUE)
b_dist_sd     <- sd(full_data_strip$avg_distance,      na.rm = TRUE)
b_height_mean <- mean(full_data_strip$avg_norm_height, na.rm = TRUE)
b_height_sd   <- sd(full_data_strip$avg_norm_height,   na.rm = TRUE)
b_hit_mean    <- mean(full_data_strip$hit_rate,        na.rm = TRUE)
b_hit_sd      <- sd(full_data_strip$hit_rate,          na.rm = TRUE)

# Initiation — from scale attributes
f_dist_mean   <- attr(decision_data_foraging_init$next_avg_distance_z,    "scaled:center")
f_dist_sd     <- attr(decision_data_foraging_init$next_avg_distance_z,    "scaled:scale")
f_height_mean <- attr(decision_data_foraging_init$next_avg_norm_height_z, "scaled:center")
f_height_sd   <- attr(decision_data_foraging_init$next_avg_norm_height_z, "scaled:scale")
f_hit_mean    <- attr(decision_data_foraging_init$next_hit_rate_z,        "scaled:center")
f_hit_sd      <- attr(decision_data_foraging_init$next_hit_rate_z,        "scaled:scale")

t_dist_mean   <- attr(decision_data_traverse_init$next_avg_distance_z,    "scaled:center")
t_dist_sd     <- attr(decision_data_traverse_init$next_avg_distance_z,    "scaled:scale")
t_height_mean <- attr(decision_data_traverse_init$next_avg_norm_height_z, "scaled:center")
t_height_sd   <- attr(decision_data_traverse_init$next_avg_norm_height_z, "scaled:scale")
t_hit_mean    <- attr(decision_data_traverse_init$next_hit_rate_z,        "scaled:center")
t_hit_sd      <- attr(decision_data_traverse_init$next_hit_rate_z,        "scaled:scale")

# ── Biting plots ───────────────────────────────────────────────────

p_bite_dist <- plot_brms_fit(
  bite_dist_brms, full_data_strip,
  x_var_z   = "avg_distance_z",
  x_var_raw = "avg_distance",
  x_mean    = b_dist_mean, x_sd = b_dist_sd,
  y_var     = "Seafloor_bites_binary",
  x_label   = "Average distance (m)",
  y_label   = "Probability of biting",
  social_var = "social_state"
)

p_bite_height <- plot_brms_fit(
  bite_height_brms, full_data_strip,
  x_var_z   = "avg_norm_height_z",
  x_var_raw = "avg_norm_height",
  x_mean    = b_height_mean, x_sd = b_height_sd,
  y_var     = "Seafloor_bites_binary",
  x_label   = "Average normalised height (m)",
  y_label   = "Probability of biting",
  social_var = "social_state"
)

p_bite_hit <- plot_brms_fit(
  bite_hitrate_brms, full_data_strip,
  x_var_z   = "hit_rate_z",
  x_var_raw = "hit_rate",
  x_mean    = b_hit_mean, x_sd = b_hit_sd,
  y_var     = "Seafloor_bites_binary",
  x_label   = "Hit rate (%)",
  y_label   = "Probability of biting",
  social_var = "social_state"
)

# ── Foraging initiation plots ──────────────────────────────────────

p_finit_dist <- plot_brms_fit(
  f_init_dist_brms, decision_data_foraging_init,
  x_var_z   = "next_avg_distance_z",
  x_var_raw = "next_avg_distance",
  x_mean    = f_dist_mean, x_sd = f_dist_sd,
  y_var     = "next_state_different",
  x_label   = "Average distance (m)",
  y_label   = "Probability of initiating foraging",
  social_var = "next_social_state"
)

p_finit_height <- plot_brms_fit(
  f_init_height_brms, decision_data_foraging_init,
  x_var_z   = "next_avg_norm_height_z",
  x_var_raw = "next_avg_norm_height",
  x_mean    = f_height_mean, x_sd = f_height_sd,
  y_var     = "next_state_different",
  x_label   = "Average normalised height (m)",
  y_label   = "Probability of initiating foraging",
  social_var = "next_social_state"
)

p_finit_hit <- plot_brms_fit(
  f_init_hitrate_brms, decision_data_foraging_init,
  x_var_z   = "next_hit_rate_z",
  x_var_raw = "next_hit_rate",
  x_mean    = f_hit_mean, x_sd = f_hit_sd,
  y_var     = "next_state_different",
  x_label   = "Hit rate (%)",
  y_label   = "Probability of initiating foraging",
  social_var = "next_social_state"
)

# ── Traversing initiation plots ────────────────────────────────────

p_tinit_dist <- plot_brms_fit(
  t_init_dist_brms, decision_data_traverse_init,
  x_var_z   = "next_avg_distance_z",
  x_var_raw = "next_avg_distance",
  x_mean    = t_dist_mean, x_sd = t_dist_sd,
  y_var     = "next_state_different",
  x_label   = "Average distance (m)",
  y_label   = "Probability of initiating traversing",
  social_var = "next_social_state"
)

p_tinit_height <- plot_brms_fit(
  t_init_height_brms, decision_data_traverse_init,
  x_var_z   = "next_avg_norm_height_z",
  x_var_raw = "next_avg_norm_height",
  x_mean    = t_height_mean, x_sd = t_height_sd,
  y_var     = "next_state_different",
  x_label   = "Average normalised height (m)",
  y_label   = "Probability of initiating traversing",
  social_var = "next_social_state"
)

p_tinit_hit <- plot_brms_fit(
  t_init_hitrate_brms, decision_data_traverse_init,
  x_var_z   = "next_hit_rate_z",
  x_var_raw = "next_hit_rate",
  x_mean    = t_hit_mean, x_sd = t_hit_sd,
  y_var     = "next_state_different",
  x_label   = "Hit rate (%)",
  y_label   = "Probability of initiating traversing",
  social_var = "next_social_state"
)

# ── Shared legend ──────────────────────────────────────────────────

# Get legend from one plot
legend_plot <- plot_brms_fit(
  bite_dist_brms, full_data_strip,
  x_var_z    = "avg_distance_z",
  x_var_raw  = "avg_distance",
  x_mean     = b_dist_mean, x_sd = b_dist_sd,
  y_var      = "Seafloor_bites_binary",
  x_label    = "Average distance (m)",
  y_label    = "Probability of biting",
  social_var = "social_state"
) + theme(legend.position = "bottom")

legend <- cowplot::get_legend(legend_plot)

# ── Combine ────────────────────────────────────────────────────────

library(cowplot)

main_plot <- plot_grid(
  p_bite_dist,   p_bite_height,   p_bite_hit,
  p_finit_dist,  p_finit_height,  p_finit_hit,
  p_tinit_dist,  p_tinit_height,  p_tinit_hit,
  ncol = 3
)

#main_plot<- plot_grid(p_bite_dist,   p_bite_height,   p_bite_hit, ncol = 3)

g <- plot_grid(main_plot, legend,
          ncol    = 1,
          rel_heights = c(1, 0.05))
g
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/brms_metric_compare_bites.png", plot = g, width = 25, height = 15)


#if time - movement
#step length
step_length_dist_current_mod <- brm(
  step_length ~ avg_distance_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = lognormal(),
  control = list(adapt_delta = 0.98)
)



step_length_height_current_mod <- brm(
  step_length ~ avg_norm_height_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = lognormal(),
  control = list(adapt_delta = 0.98)
)

summary(step_length_height_current_mod)
step_length_hitrate_current_mod <- brm(
  step_length ~ hit_rate_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = lognormal(),
  control = list(adapt_delta = 0.98)
)

#turning angle
hist(full_movement_data$turn_angle_rad_3d)
class(full_movement_data$turn_angle_rad_3d)
full_movement_data_TA <- full_movement_data %>%
  tidyr::drop_na(turn_angle_rad_3d)

turn_angle_rad_3d_dist_current_mod <- brm(
  turn_angle_rad_3d ~ avg_distance_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data_TA,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = lognormal(), #may change to beta regression as bounded under 3.14
  control = list(adapt_delta = 0.98)
)

turn_angle_rad_3d_height_current_mod <- brm(
  turn_angle_rad_3d ~ avg_norm_height_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data_TA,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = lognormal(), #may change to beta regression as bounded under 3.14
  control = list(adapt_delta = 0.98)
)

turn_angle_rad_3d_hitrate_current_mod <- brm(
  turn_angle_rad_3d ~ hit_rate_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data_TA,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = lognormal(), #may change to beta regression as bounded under 3.14
  control = list(adapt_delta = 0.98)
)

#straigntness window
full_movement_data_St <- full_movement_data %>%
  tidyr::drop_na(straightness3D)
straightness3D_dist_current_mod <- brm(
  straightness3D ~ avg_distance_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data_St,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(gamma(2,0.1), class = "phi")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = Beta(),
  control = list(adapt_delta = 0.98)
)

summary(straightness3D_dist_current_mod)

straightness3D_height_current_mod <- brm(
  straightness3D ~ avg_norm_height_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data_St,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(gamma(2,0.1), class = "phi")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = Beta(),
  control = list(adapt_delta = 0.98)
)

summary(straightness3D_height_current_mod)

straightness3D_hitrate_current_mod <- brm(
  straightness3D ~ hit_rate_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = full_movement_data_St,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    prior(normal(0, 1), class = "Intercept"),  
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(gamma(2,0.1), class = "phi")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = Beta(),
  control = list(adapt_delta = 0.98)
)

# ── Scale attributes ───────────────────────────────────────────────

m_dist_mean   <- 0.3768193
m_dist_sd     <- 0.3010322
m_height_mean <- 0.1503749
m_height_sd   <- 0.1252082
m_hit_mean    <- 52.91566
m_hit_sd      <- 10.34143

# ── Generic plot function for continuous outcomes ──────────────────

plot_brms_continuous <- function(model, data,
                                 x_var_z,
                                 x_var_raw,
                                 x_mean, x_sd,
                                 y_var,
                                 x_label, y_label,
                                 social_var = "social_state",
                                 ndraws = 500,
                                 n_bins = 20) {
  
  # ── Binned raw data ──────────────────────────────────────────────
  binned <- data %>%
    filter(!is.na(.data[[x_var_raw]]) & !is.na(.data[[y_var]])) %>%
    mutate(x_bin = cut(.data[[x_var_raw]], breaks = n_bins)) %>%
    group_by(x_bin, subject, .data[[social_var]]) %>%
    summarise(
      x_mid  = mean(.data[[x_var_raw]], na.rm = TRUE),
      y_mean = mean(.data[[y_var]],     na.rm = TRUE),
      n      = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(x_mid))
  
  # ── Per-combination prediction grid ─────────────────────────────
  x_ranges <- data %>%
    filter(!is.na(.data[[x_var_raw]])) %>%
    group_by(subject, .data[[social_var]]) %>%
    summarise(
      x_min = min(.data[[x_var_raw]], na.rm = TRUE),
      x_max = max(.data[[x_var_raw]], na.rm = TRUE),
      .groups = "drop"
    )
  
  pred_list <- vector("list", nrow(x_ranges))
  for (i in seq_len(nrow(x_ranges))) {
    pred_list[[i]] <- data.frame(
      subject          = x_ranges$subject[i],
      social           = x_ranges[[social_var]][i],
      x_raw            = seq(x_ranges$x_min[i],
                             x_ranges$x_max[i],
                             length.out = 100),
      stringsAsFactors = FALSE
    )
  }
  
  pred_grid <- bind_rows(pred_list)
  names(pred_grid)[names(pred_grid) == "social"] <- social_var
  pred_grid[[x_var_z]] <- (pred_grid$x_raw - x_mean) / x_sd
  pred_grid$Trial_ID    <- NA
  pred_grid$site        <- NA
  pred_grid$folder_date <- NA
  
  # ── Posterior predictions ────────────────────────────────────────
  pred_draws <- pred_grid %>%
    add_epred_draws(model, re_formula = NA, ndraws = ndraws)
  
  pred_summary <- pred_draws %>%
    group_by(x_raw, subject, .data[[social_var]]) %>%
    summarise(
      median = median(.epred),
      lower  = quantile(.epred, 0.025),
      upper  = quantile(.epred, 0.975),
      .groups = "drop"
    )
  
  # ── Plot ─────────────────────────────────────────────────────────
  ggplot() +
    geom_point(
      data  = binned,
      aes(x     = x_mid,
          y     = y_mean,
          color = .data[[social_var]],
          size  = n),
      alpha = 0.5
    ) +
    geom_ribbon(
      data = pred_summary,
      aes(x    = x_raw,
          ymin = lower,
          ymax = upper,
          fill = .data[[social_var]]),
      alpha = 0.2
    ) +
    geom_line(
      data = pred_summary,
      aes(x     = x_raw,
          y     = median,
          color = .data[[social_var]]),
      linewidth = 1
    ) +
    scale_size_continuous(name = "N", range = c(1, 6)) +
    facet_wrap(~ subject) +
    scale_color_manual(
      values = c("Asocial" = "#F8766D",
                 "Social"  = "#00BFC4")
    ) +
    scale_fill_manual(
      values = c("Asocial" = "#F8766D",
                 "Social"  = "#00BFC4")
    ) +
    labs(x     = x_label,
         y     = y_label,
         color = "Social state",
         fill  = "Social state",
         size  = "N") +
    theme_classic(base_size = 18) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          legend.position  = "none")
}

# ── Step length plots ──────────────────────────────────────────────

p_sl_dist <- plot_brms_continuous(
  step_length_dist_current_mod, full_movement_data,
  x_var_z   = "avg_distance_z",
  x_var_raw = "avg_distance",
  x_mean    = m_dist_mean, x_sd = m_dist_sd,
  y_var     = "step_length",
  x_label   = "Average distance (m)",
  y_label   = "Step length (m)"
)

p_sl_height <- plot_brms_continuous(
  step_length_height_current_mod, full_movement_data,
  x_var_z   = "avg_norm_height_z",
  x_var_raw = "avg_norm_height",
  x_mean    = m_height_mean, x_sd = m_height_sd,
  y_var     = "step_length",
  x_label   = "Average normalised height (m)",
  y_label   = "Step length (m)"
)

p_sl_hit <- plot_brms_continuous(
  step_length_hitrate_current_mod, full_movement_data,
  x_var_z   = "hit_rate_z",
  x_var_raw = "hit_rate",
  x_mean    = m_hit_mean, x_sd = m_hit_sd,
  y_var     = "step_length",
  x_label   = "Hit rate (%)",
  y_label   = "Step length (m)"
)

# ── Turning angle plots ────────────────────────────────────────────

p_ta_dist <- plot_brms_continuous(
  turn_angle_rad_3d_dist_current_mod, full_movement_data_TA,
  x_var_z   = "avg_distance_z",
  x_var_raw = "avg_distance",
  x_mean    = m_dist_mean, x_sd = m_dist_sd,
  y_var     = "turn_angle_rad_3d",
  x_label   = "Average distance (m)",
  y_label   = "Turning angle (rad)"
)

p_ta_height <- plot_brms_continuous(
  turn_angle_rad_3d_height_current_mod, full_movement_data_TA,
  x_var_z   = "avg_norm_height_z",
  x_var_raw = "avg_norm_height",
  x_mean    = m_height_mean, x_sd = m_height_sd,
  y_var     = "turn_angle_rad_3d",
  x_label   = "Average normalised height (m)",
  y_label   = "Turning angle (rad)"
)

p_ta_hit <- plot_brms_continuous(
  turn_angle_rad_3d_hitrate_current_mod, full_movement_data_TA,
  x_var_z   = "hit_rate_z",
  x_var_raw = "hit_rate",
  x_mean    = m_hit_mean, x_sd = m_hit_sd,
  y_var     = "turn_angle_rad_3d",
  x_label   = "Hit rate (%)",
  y_label   = "Turning angle (rad)"
)

# ── Straightness plots ─────────────────────────────────────────────

p_st_dist <- plot_brms_continuous(
  straightness3D_dist_current_mod, full_movement_data_St,
  x_var_z   = "avg_distance_z",
  x_var_raw = "avg_distance",
  x_mean    = m_dist_mean, x_sd = m_dist_sd,
  y_var     = "straightness3D",
  x_label   = "Average distance (m)",
  y_label   = "Straightness"
)

p_st_height <- plot_brms_continuous(
  straightness3D_height_current_mod, full_movement_data_St,
  x_var_z   = "avg_norm_height_z",
  x_var_raw = "avg_norm_height",
  x_mean    = m_height_mean, x_sd = m_height_sd,
  y_var     = "straightness3D",
  x_label   = "Average normalised height (m)",
  y_label   = "Straightness"
)

p_st_hit <- plot_brms_continuous(
  straightness3D_hitrate_current_mod, full_movement_data_St,
  x_var_z   = "hit_rate_z",
  x_var_raw = "hit_rate",
  x_mean    = m_hit_mean, x_sd = m_hit_sd,
  y_var     = "straightness3D",
  x_label   = "Hit rate (%)",
  y_label   = "Straightness"
)

# ── Shared legend ──────────────────────────────────────────────────

legend_plot <- plot_brms_continuous(
  step_length_dist_current_mod, full_movement_data,
  x_var_z   = "avg_distance_z",
  x_var_raw = "avg_distance",
  x_mean    = m_dist_mean, x_sd = m_dist_sd,
  y_var     = "step_length",
  x_label   = "Average distance (m)",
  y_label   = "Step length (m)"
) + theme(legend.position = "bottom")

legend <- cowplot::get_legend(legend_plot)

# ── Combine ────────────────────────────────────────────────────────

main_plot <- cowplot::plot_grid(
  p_sl_dist,   p_sl_height,   p_sl_hit,
  p_ta_dist,   p_ta_height,   p_ta_hit,
  p_st_dist,   p_st_height,   p_st_hit,
  ncol = 3
)

g <- cowplot::plot_grid(
  main_plot, legend,
  ncol        = 1,
  rel_heights = c(1, 0.05)
)
g
ggsave(
  "/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/movement_metrics.png",
  plot   = g,
  width  = 22,
  height = 14
)


##log scale for behavior

# ── Biting models ──────────────────────────────────────────────────

bite_log_dist_brms <- brm(
  Seafloor_bites_binary ~ log_avg_distance_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = full_data_strip,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),    class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),  class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "bite_log_dist_brms"  # saves model so you don't refit
)

summary(bite_dist_brms)

bite_log_height_brms <- brm(
  Seafloor_bites_binary ~ log_avg_norm_height_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = full_data_strip,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "bite_log_height_brms"
)

summary(bite_height_brms)

# ── Foraging initiation models ─────────────────────────────────────

f_init_log_dist_brms <- brm(
  next_state_different ~ log_next_avg_distance_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_foraging_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "finit_log_dist_brms"
)


f_init_log_height_brms <- brm(
  next_state_different ~ log_next_avg_norm_height_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_foraging_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "finit_log_height_brms"
)
summary(f_init_height_brms)


# ── Traversing initiation models ─────────────────────────────────────
t_init_log_dist_brms <- brm(
  next_state_different ~ log_next_avg_distance_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_traverse_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "tinit_log_dist_brms"
)

t_init_log_height_brms <- brm(
  next_state_different ~ log_next_avg_norm_height_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data    = decision_data_traverse_init,
  family  = brms::bernoulli(),
  chains  = 4, cores = 4,
  iter    = 4000, warmup = 2000,
  prior   = c(
    prior(normal(0, 1),      class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1),    class = "sd")
  ),
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.99),
  file    = "tinit_log_height_brms"
)
summary(t_init_height_brms)


###############
# ── Scale attributes ───────────────────────────────────────────────

# Biting log scale
b_log_dist_mean   <- -1.211119
b_log_dist_sd     <-  0.6594093
b_log_height_mean <- -2.153199
b_log_height_sd   <-  0.6965962

# Foraging initiation log scale
f_log_dist_mean   <- -1.210903
f_log_dist_sd     <-  0.6593948
f_log_height_mean <- -2.153191
f_log_height_sd   <-  0.6967062

# Traversing initiation — check these
t_log_dist_mean   <- attr(decision_data_traverse_init$log_next_avg_distance_z,    "scaled:center")
t_log_dist_sd     <- attr(decision_data_traverse_init$log_next_avg_distance_z,    "scaled:scale")
t_log_height_mean <- attr(decision_data_traverse_init$log_next_avg_norm_height_z, "scaled:center")
t_log_height_sd   <- attr(decision_data_traverse_init$log_next_avg_norm_height_z, "scaled:scale")

# ── Plot function for log-scale models ────────────────────────────
# x axis shows raw values, model was fit on log scale

plot_brms_log <- function(model, data,
                          x_var_z,        # log scaled variable in model
                          x_var_raw,      # raw (not log) variable for x axis
                          x_log_mean,     # mean of log variable
                          x_log_sd,       # sd of log variable
                          y_var,
                          x_label, y_label,
                          social_var,
                          ndraws = 500,
                          n_bins = 20) {
  
  # ── Raw counts at 0 and 1 ────────────────────────────────────────
  raw_counts <- data %>%
    filter(!is.na(.data[[x_var_raw]])) %>%
    mutate(x_bin = cut(.data[[x_var_raw]], breaks = n_bins)) %>%
    group_by(x_bin, subject, .data[[social_var]], .data[[y_var]]) %>%
    summarise(
      x_mid = mean(.data[[x_var_raw]], na.rm = TRUE),
      n     = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(x_mid)) %>%
    mutate(
      y_offset = case_when(
        .data[[social_var]] == "Asocial" & .data[[y_var]] == 0 ~ -0.02,
        .data[[social_var]] == "Asocial" & .data[[y_var]] == 1 ~  0.98,
        .data[[social_var]] == "Social"  & .data[[y_var]] == 0 ~  0.02,
        .data[[social_var]] == "Social"  & .data[[y_var]] == 1 ~  1.02
      )
    )
  
  # ── Per-combination prediction grid ──────────────────────────────
  # Build on raw scale, then log-transform and scale for prediction
  x_ranges <- data %>%
    filter(!is.na(.data[[x_var_raw]]) & .data[[x_var_raw]] > 0) %>%
    group_by(subject, .data[[social_var]]) %>%
    summarise(
      x_min = min(.data[[x_var_raw]], na.rm = TRUE),
      x_max = max(.data[[x_var_raw]], na.rm = TRUE),
      .groups = "drop"
    )
  
  pred_list <- vector("list", nrow(x_ranges))
  for (i in seq_len(nrow(x_ranges))) {
    pred_list[[i]] <- data.frame(
      subject          = x_ranges$subject[i],
      social           = x_ranges[[social_var]][i],
      x_raw            = seq(x_ranges$x_min[i],
                             x_ranges$x_max[i],
                             length.out = 100),
      stringsAsFactors = FALSE
    )
  }
  
  pred_grid <- bind_rows(pred_list)
  names(pred_grid)[names(pred_grid) == "social"] <- social_var
  
  # Log transform then scale — matches how model was fitted
  pred_grid[[x_var_z]]  <- (log(pred_grid$x_raw) - x_log_mean) / x_log_sd
  pred_grid$Trial_ID    <- NA
  pred_grid$site        <- NA
  pred_grid$folder_date <- NA
  
  if ("time_in_state" %in% colnames(data)) {
    pred_grid$time_in_state <- median(data$time_in_state, na.rm = TRUE)
  }
  
  # ── Posterior predictions ────────────────────────────────────────
  pred_draws <- pred_grid %>%
    add_epred_draws(model, re_formula = NA, ndraws = ndraws)
  
  pred_summary <- pred_draws %>%
    group_by(x_raw, subject, .data[[social_var]]) %>%
    summarise(
      median = median(.epred),
      lower  = quantile(.epred, 0.025),
      upper  = quantile(.epred, 0.975),
      .groups = "drop"
    )
  
  # ── Plot ─────────────────────────────────────────────────────────
  ggplot() +
    geom_point(
      data  = raw_counts,
      aes(x     = x_mid,
          y     = y_offset,
          color = .data[[social_var]],
          size  = n),
      alpha = 0.4
    ) +
    geom_ribbon(
      data  = pred_summary,
      aes(x    = x_raw,
          ymin = lower,
          ymax = upper,
          fill = .data[[social_var]]),
      alpha = 0.2
    ) +
    geom_line(
      data = pred_summary,
      aes(x     = x_raw,
          y     = median,
          color = .data[[social_var]]),
      linewidth = 1
    ) +
    scale_size_continuous(name = "N", range = c(0.5, 6)) +
    facet_wrap(~ subject) +
    scale_color_manual(
      values = c("Asocial" = "#F8766D",
                 "Social"  = "#00BFC4")
    ) +
    scale_fill_manual(
      values = c("Asocial" = "#F8766D",
                 "Social"  = "#00BFC4")
    ) +
    labs(x     = x_label,
         y     = y_label,
         color = "Social state",
         fill  = "Social state",
         size  = "N") +
    theme_classic(base_size = 18) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          legend.position  = "none")
}

# ── Biting log plots ───────────────────────────────────────────────

p_bite_log_dist <- plot_brms_log(
  bite_log_dist_brms, full_data_strip,
  x_var_z       = "log_avg_distance_z",
  x_var_raw     = "avg_distance",
  x_log_mean    = b_log_dist_mean,
  x_log_sd      = b_log_dist_sd,
  y_var         = "Seafloor_bites_binary",
  x_label       = "Average distance (m)",
  y_label       = "Probability of biting",
  social_var    = "social_state"
)

p_bite_log_height <- plot_brms_log(
  bite_log_height_brms, full_data_strip,
  x_var_z       = "log_avg_norm_height_z",
  x_var_raw     = "avg_norm_height",
  x_log_mean    = b_log_height_mean,
  x_log_sd      = b_log_height_sd,
  y_var         = "Seafloor_bites_binary",
  x_label       = "Average normalised height (m)",
  y_label       = "Probability of biting",
  social_var    = "social_state"
)

# ── Foraging initiation log plots ──────────────────────────────────

p_finit_log_dist <- plot_brms_log(
  f_init_log_dist_brms, decision_data_foraging_init,
  x_var_z       = "log_next_avg_distance_z",
  x_var_raw     = "next_avg_distance",
  x_log_mean    = f_log_dist_mean,
  x_log_sd      = f_log_dist_sd,
  y_var         = "next_state_different",
  x_label       = "Average distance (m)",
  y_label       = "Probability of initiating foraging",
  social_var    = "next_social_state"
)

p_finit_log_height <- plot_brms_log(
  f_init_log_height_brms, decision_data_foraging_init,
  x_var_z       = "log_next_avg_norm_height_z",
  x_var_raw     = "next_avg_norm_height",
  x_log_mean    = f_log_height_mean,
  x_log_sd      = f_log_height_sd,
  y_var         = "next_state_different",
  x_label       = "Average normalised height (m)",
  y_label       = "Probability of initiating foraging",
  social_var    = "next_social_state"
)

# ── Traversing initiation log plots ───────────────────────────────

p_tinit_log_dist <- plot_brms_log(
  t_init_log_dist_brms, decision_data_traverse_init,
  x_var_z       = "log_next_avg_distance_z",
  x_var_raw     = "next_avg_distance",
  x_log_mean    = t_log_dist_mean,
  x_log_sd      = t_log_dist_sd,
  y_var         = "next_state_different",
  x_label       = "Average distance (m)",
  y_label       = "Probability of initiating traversing",
  social_var    = "next_social_state"
)

p_tinit_log_height <- plot_brms_log(
  t_init_log_height_brms, decision_data_traverse_init,
  x_var_z       = "log_next_avg_norm_height_z",
  x_var_raw     = "next_avg_norm_height",
  x_log_mean    = t_log_height_mean,
  x_log_sd      = t_log_height_sd,
  y_var         = "next_state_different",
  x_label       = "Average normalised height (m)",
  y_label       = "Probability of initiating traversing",
  social_var    = "next_social_state"
)

# ── Combine 2 columns × 3 rows ────────────────────────────────────

main_plot <- cowplot::plot_grid(
  p_bite_log_dist,   p_bite_log_height,
  p_finit_log_dist,  p_finit_log_height,
  p_tinit_log_dist,  p_tinit_log_height,
  ncol = 2
)

legend_plot <- plot_brms_log(
  bite_log_dist_brms, full_data_strip,
  x_var_z    = "log_avg_distance_z",
  x_var_raw  = "avg_distance",
  x_log_mean = b_log_dist_mean,
  x_log_sd   = b_log_dist_sd,
  y_var      = "Seafloor_bites_binary",
  x_label    = "Average distance (m)",
  y_label    = "Probability of biting",
  social_var = "social_state"
) + theme(legend.position = "bottom")

legend <- cowplot::get_legend(legend_plot)

g <- cowplot::plot_grid(
  main_plot, legend,
  ncol        = 1,
  rel_heights = c(1, 0.05)
)

g

ggsave(
  "/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/brms_log_behavior.png",
  plot   = g,
  width  = 16,
  height = 18
)
