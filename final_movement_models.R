
#functions
# Helper function to compute straightness over a window
sliding_straightness <- function(x, y, z, width = 5) {
  n <- length(x)
  half <- floor(width / 2)
  result <- rep(NA_real_, n)
  
  for (i in (half + 1):(n - half)) {
    idx <- (i - half):(i + half)
    result[i] <- calculate_straightness(x[idx], y[idx], z[idx])
  }
  result
}
calculate_straightness <- function(x, y, z) {
  # Straightness index: net displacement / total path length
  # Ranges from 0 (very tortuous) to 1 (perfectly straight)
  
  n <- length(x)
  if (n < 2 || all(is.na(x)) || all(is.na(y))) return(NA_real_)
  
  # Net displacement (straight-line distance from first to last point)
  net_displacement <- sqrt((x[n]-x[1])^2 + (y[n]-y[1])^2 + (z[n]-z[1])^2)
  
  # Total path length (sum of step lengths)
  step_lengths <- sqrt(diff(x)^2 + diff(y)^2 + diff(z)^2)
  total_path <- sum(step_lengths, na.rm = TRUE)
  
  if (total_path == 0) return(NA_real_)
  
  net_displacement / total_path
}

#get movement in
full_data$x <- as.numeric(full_data$x)
full_data$y <- as.numeric(full_data$y)
full_data$z <- as.numeric(full_data$z)

full_movement_data <- full_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  filter(!is.na(x)) %>%
  mutate(
    dX = lead(x) - x,
    dY = lead(y) - y,
    dZ = lead(z) - z,
    dt = lead(sec) - sec,
    step_length = sqrt(dX^2 + dY^2 + dZ^2), # 3D step length
    next_log_ehc = lead(log_ehc),
    next_ehc = lead(ehc),
    next_social_state = lead(social_state),
    next_activity_state = lead(activity_state),
    dot_product = dX * lag(dX) + dY * lag(dY) + dZ * lag(dZ),
    mag_current = step_length,
    mag_prev = lag(step_length),
    turn_angle_rad_3d = acos( # dot product of consecutive step vectors divided by product of their magnitudes
      pmin(pmax(dot_product / (mag_current * mag_prev), -1), 1)  # clamp to [-1,1] to avoid NaN from floating point
    ),
    bearing = atan2(dY, dX), # 2D Bearing of each step (in radians, atan2 handles all quadrants). Has left and right turn?
    turn_angle_rad = bearing - lag(bearing), # Turning angle = change in bearing between successive steps
    turn_angle_rad = atan2(sin(turn_angle_rad), cos(turn_angle_rad)), # Wrap to (-pi, pi] so angles don't exceed + or - 180
    elevation_angle = asin( #(pitch): angle above/below horizontal plane
      pmin(pmax(dZ / step_length, -1), 1)  # clamp to [-1,1]
    ),
    speed = step_length / dt, # Speed
    straightness3D = sliding_straightness(x, y, z, width = 5)# straightness indices over 5, 10, and 15 s window     # 2D straightness index (5 second window)
  )  %>%
  ungroup() %>%
  filter(dt == 1 & !is.na(next_log_ehc)) %>%
  mutate(
    next_log_ehc_z = scale(next_log_ehc),
    next_ehc_z = scale(next_ehc),
    next_social_state = factor(next_social_state)
  )



full_movement_data %>%
  ggplot(aes(x = next_social_state, y = straightness3D)) +
  geom_boxplot() +
  geom_smooth()


#Step length with current
step_length_current_mod <- brm(
  step_length ~ log_ehc_z*social_state*subject +
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
saveRDS(step_length_current_mod, file = "step_length_current_mod.rds")

pp_check(step_length_current_mod) +
  ggplot2::scale_x_log10()

summary(step_length_current_mod)

saveRDS(step_length_current_mod, file = "step_length_current_mod.rds")

#fig
pred_newd <- full_movement_data |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

observed_combinations <- full_movement_data |>
  dplyr::distinct(subject, social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_state")
  )

pmu <- posterior_epred(
  step_length_current_mod,
  newdata = pred_newd,
  re_formula = NA
)

preds <- cbind(
  pred_newd,
  mnmu = apply(pmu, 2, mean),
  t(apply(pmu, 2, hpdi, prob = 0.95))
)

names(preds)[names(preds) == "1"] <- "mulo95"
names(preds)[names(preds) == "2"] <- "muhi95"

orig_mean <- attr(full_movement_data$log_ehc_z, "scaled:center")
orig_sd <- attr(full_movement_data$log_ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

Plot
class()
step_length_current_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = full_movement_data,
    aes(
      x = log_ehc,
      y = step_length,
      color = social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = log_ehc,
      ymin = mulo95,
      ymax = muhi95,
      fill = social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = log_ehc,
      y = mnmu,
      color = social_state
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Step length"
  ) +
  theme_classic(base_size = 24)

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/step_length_current_fig.png", plot = step_length_current_fig, width = 10, height = 5)



#Step length with lead
hist(log(full_movement_data$step_length))
hist(full_movement_data$next_ehc_z)

class(full_movement_data$next_social_state)
step_length_lead_mod <- brm(
  step_length ~ next_log_ehc_z*next_social_state*subject +
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
  control = list(adapt_delta = 0.99)
)
saveRDS(step_length_lead_mod, file = "step_length_lead_mod.rds")
pp_check(step_length_lead_mod, ndraws = 100)
summary(step_length_current_mod)


#fig
pred_newd <- full_movement_data |>
  dplyr::group_by(subject, next_social_state) |>
  dplyr::summarise(
    next_log_ehc_z = list(seq(min(next_log_ehc_z), max(next_log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(next_log_ehc_z)

observed_combinations <- full_movement_data |>
  dplyr::distinct(subject, next_social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "next_social_state")
  )

pmu <- posterior_epred(
  step_length_lead_mod,
  newdata = pred_newd,
  re_formula = NA
)

preds <- cbind(
  pred_newd,
  mnmu = apply(pmu, 2, mean),
  t(apply(pmu, 2, hpdi, prob = 0.95))
)

names(preds)[names(preds) == "1"] <- "mulo95"
names(preds)[names(preds) == "2"] <- "muhi95"

orig_mean <- attr(full_movement_data$next_log_ehc_z, "scaled:center")
orig_sd <- attr(full_movement_data$next_log_ehc_z, "scaled:scale")

preds$next_log_ehc <- preds$next_log_ehc_z * orig_sd + orig_mean


step_length_lead_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = full_movement_data,
    aes(
      x = next_log_ehc,
      y = step_length,
      color = next_social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = next_log_ehc,
      ymin = mulo95,
      ymax = muhi95,
      fill = next_social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = next_log_ehc,
      y = mnmu,
      color = next_social_state
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Step length"
  ) +
  theme_classic(base_size = 24)

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/step_length_lead_fig.png", plot = step_length_lead_fig, width = 10, height = 5)

step_length_activity_state_fig <- full_movement_data %>%
  ggplot() +
  geom_histogram(aes(x = log(step_length), fill = as.factor(activity_state)), position = "identity", alpha = 0.5) +
  facet_wrap(next_social_state~subject) +
  theme_minimal()
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/step_length_activity_state_fig.png", plot = step_length_activity_state_fig, width = 10, height = 5)

step_length_bite_fig <- full_movement_data %>%
  ggplot() +
  geom_histogram(aes(x = log(step_length), fill = as.factor(Seafloor_bites_binary)), position = "identity", alpha = 0.5) +
  facet_wrap(next_social_state~subject) +
  theme_minimal()
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/step_length_bites_fig.png", plot = step_length_bite_fig, width = 10, height = 5)

step_length_activity_state_density_fig <- full_movement_data %>%
  ggplot() +
  geom_density(aes(x = log(step_length), fill = as.factor(activity_state)), alpha = 0.5) +
  facet_wrap(next_social_state~subject) +
  theme_minimal()
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/step_length_activity_state_fig.png", plot = step_length_activity_state_density_fig, width = 10, height = 5)

step_length_bite_density_fig <- full_movement_data %>%
  ggplot() +
  geom_density(aes(x = log(step_length), fill = as.factor(Seafloor_bites_binary)), alpha = 0.5) +
  facet_wrap(next_social_state~subject) +
  theme_minimal()
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/step_length_bites_fig.png", plot = step_length_bite_fig, width = 10, height = 5)


#turn angle rad 3d with current
hist(full_movement_data$turn_angle_rad_3d)
class(full_movement_data$turn_angle_rad_3d)
full_movement_data_TA <- full_movement_data %>%
  tidyr::drop_na(turn_angle_rad_3d)

turn_angle_rad_3d_current_mod <- brm(
  turn_angle_rad_3d ~ log_ehc_z*social_state*subject +
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
summary(turn_angle_rad_3d_current_mod)
pp_check(turn_angle_rad_3d_current_mod, ndraws = 100) +
  ggplot2::scale_x_log10()
saveRDS(turn_angle_rad_3d_current_mod, file = "turn_angle_rad_3d_current_mod.rds")

#fig
pred_newd <- full_movement_data |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

observed_combinations <- full_movement_data |>
  dplyr::distinct(subject, social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_state")
  )

pmu <- posterior_epred(
  turn_angle_rad_3d_current_mod,
  newdata = pred_newd,
  re_formula = NA
)

preds <- cbind(
  pred_newd,
  mnmu = apply(pmu, 2, mean),
  t(apply(pmu, 2, hpdi, prob = 0.95))
)

names(preds)[names(preds) == "1"] <- "mulo95"
names(preds)[names(preds) == "2"] <- "muhi95"

orig_mean <- attr(full_movement_data$log_ehc_z, "scaled:center")
orig_sd <- attr(full_movement_data$log_ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean


turn_angle_rad_3d_current_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = full_movement_data,
    aes(
      x = log_ehc,
      y = turn_angle_rad_3d,
      color = social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = log_ehc,
      ymin = mulo95,
      ymax = muhi95,
      fill = social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = log_ehc,
      y = mnmu,
      color = social_state
    ),
    linewidth = 1
  ) +
  coord_cartesian(ylim = c(0, max(full_movement_data$turn_angle_rad_3d, na.rm=TRUE))) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Turning angle (radians)"
  ) +
  theme_classic(base_size = 24)
turn_angle_rad_3d_current_fig
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/turn_angle_rad_3d_current_fig.png", plot = turn_angle_rad_3d_current_fig, width = 10, height = 5)


#turn angle rad 3d with lead

#turn angle rad 2d with current
#turn angle rad 2d with lead

#elevation angle with current
#elevation angle with lead


#straightness with sliding window 5
hist(full_movement_data$straightness3D)
min(full_movement_data$straightness3D, na.rm=TRUE)

full_movement_data_St <- full_movement_data %>%
  tidyr::drop_na(straightness3D)
straightness3D_current_mod <- brm(
  straightness3D ~ log_ehc_z*social_state*subject +
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

pp_check(straightness3D_current_mod)

summary(straightness3D_current_mod)
saveRDS(straightness3D_current_mod, file = "straightness3D_current_mod.rds")

#fig
pred_newd <- full_movement_data |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

observed_combinations <- full_movement_data |>
  dplyr::distinct(subject, social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_state")
  )

pmu <- posterior_epred(
  straightness3D_current_mod,
  newdata = pred_newd,
  re_formula = NA
)

preds <- cbind(
  pred_newd,
  mnmu = apply(pmu, 2, mean),
  t(apply(pmu, 2, hpdi, prob = 0.95))
)

names(preds)[names(preds) == "1"] <- "mulo95"
names(preds)[names(preds) == "2"] <- "muhi95"

orig_mean <- attr(full_movement_data$log_ehc_z, "scaled:center")
orig_sd <- attr(full_movement_data$log_ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean


straightness3D_current_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = full_movement_data,
    aes(
      x = log_ehc,
      y = straightness3D,
      color = social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = log_ehc,
      ymin = mulo95,
      ymax = muhi95,
      fill = social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = log_ehc,
      y = mnmu,
      color = social_state
    ),
    linewidth = 1
  ) +
  #coord_cartesian(ylim = c(0, max(full_movement_data$straightness3D, na.rm=TRUE))) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Straightness 3D (5s)"
  ) +
  theme_classic(base_size = 24)
straightness3D_current_fig
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/straightness3D_current_fig.png", plot = straightness3D_current_fig, width = 10, height = 5)





