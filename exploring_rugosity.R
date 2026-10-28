full_data$local_rugosity <- as.numeric(full_data$local_rugosity)
full_data$relief <- as.numeric(full_data$relief)

hist(full_data$local_rugosity)

#prelim bites
full_data %>%
  ggplot(aes(x=local_rugosity, y = Seafloor_bites_binary)) +
  geom_point() +
  geom_smooth(method.args = list(family = "binomial")) +
  facet_wrap(~subject)

#prelim foraging init
decision_data_foraging_init %>%
  ggplot(aes(x = next_local_rugosity, y = next_state_different)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~subject)

#foraging rugosity mod
traversing_initiation_relief_model <- brm(
  formula = next_state_different ~
    next_relief_z*subject*next_social_state + 
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = decision_data_traverse_init,
  
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  
  save_pars = save_pars(all = TRUE),
  
  #sample_prior = "only",
  
  family = brms::bernoulli(),
  
  control = list(
    adapt_delta = 0.99
  )
)

summary(foraging_initiation_relief_model)
#fig
pred_newd <- decision_data_foraging_init |>
  dplyr::group_by(subject, next_social_state) |>
  dplyr::summarise(
    next_relief_z = list(seq(min(next_relief_z), max(next_relief_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(next_relief_z)

observed_combinations <- decision_data_foraging_init |>
  dplyr::distinct(subject, next_social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "next_social_state")
  )

pmu <- posterior_epred(
  foraging_initiation_relief_model,
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

orig_mean <- attr(decision_data_foraging_init$next_relief_z, "scaled:center")
orig_sd <- attr(decision_data_foraging_init$next_relief_z, "scaled:scale")

preds$next_relief <- preds$next_relief_z * orig_sd + orig_mean

foraging_initiation_rugosity_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = decision_data_foraging_init,
    aes(
      x = next_relief,
      y = next_state_different,
      color = next_social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = next_relief,
      ymin = mulo95,
      ymax = muhi95,
      fill = next_social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = next_relief,
      y = mnmu,
      color = next_social_state
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Next local relief (m)",
    y = "Probability of initiating foraging"
  ) +
  theme_classic(base_size = 24)

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/foraging_initiation_rugosity_fig.png", plot = foraging_initiation_rugosity_fig, width = 10, height = 5)

#traversing rugosity mod
traversing_initiation_rugosity_model <- brm(
  formula = next_state_different ~
    next_local_rugosity_z*subject*next_social_state + 
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = decision_data_traverse_init,
  
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  
  save_pars = save_pars(all = TRUE),
  
  #sample_prior = "only",
  
  family = brms::bernoulli(),
  
  control = list(
    adapt_delta = 0.99
  )
)

summary(traversing_initiation_rugosity_model)
#fig
pred_newd <- decision_data_traverse_init |>
  dplyr::group_by(subject, next_social_state) |>
  dplyr::summarise(
    next_local_rugosity_z = list(seq(min(next_local_rugosity_z), max(next_local_rugosity_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(next_local_rugosity_z)

observed_combinations <- decision_data_traverse_init |>
  dplyr::distinct(subject, next_social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "next_social_state")
  )

pmu <- posterior_epred(
  traversing_initiation_rugosity_model,
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

orig_mean <- attr(decision_data_traverse_init$next_local_rugosity_z, "scaled:center")
orig_sd <- attr(decision_data_traverse_init$next_local_rugosity_z, "scaled:scale")

preds$next_local_rugosity <- preds$next_local_rugosity_z * orig_sd + orig_mean

foraging_initiation_rugosity_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = decision_data_traverse_init,
    aes(
      x = next_local_rugosity,
      y = next_state_different,
      color = next_social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = next_local_rugosity,
      ymin = mulo95,
      ymax = muhi95,
      fill = next_social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = next_local_rugosity,
      y = mnmu,
      color = next_social_state
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Next local rugosity",
    y = "Probability of initiating traversing"
  ) +
  theme_classic(base_size = 24)

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/foraging_initiation_rugosity_fig.png", plot = foraging_initiation_rugosity_fig, width = 10, height = 5)

##seafloor bite model
full_data_strip$local_rugosity_z <- scale(full_data_strip$local_rugosity)
seafloor_bites_rugosity_model <- brm(
  formula = Seafloor_bites_binary ~
    local_rugosity_z*subject*social_state +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = full_data_strip,
  
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  
  save_pars = save_pars(all = TRUE),
  
  #sample_prior = "only",
  
  family = brms::bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)
