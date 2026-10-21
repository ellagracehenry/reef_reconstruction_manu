Create the prediction grid

pred_newd <- full_movement_data |>
  dplyr::group_by(subject, next_social_state, Seafloor_bites_binary) |>
  dplyr::summarise(
    next_log_ehc_z = list(seq(min(next_log_ehc_z), max(next_log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(next_log_ehc_z)

Keep only subject × social_state combinations that actually occur in the data

observed_combinations <- full_movement_data |>
  dplyr::distinct(subject, next_social_state, Seafloor_bites_binary)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "next_social_state", "Seafloor_bites_binary")
  )


Generate predictions

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

Undo standardisation

orig_mean <- attr(full_movement_data$next_log_ehc_z, "scaled:center")
orig_sd <- attr(full_movement_data$next_log_ehc_z, "scaled:scale")

preds$next_log_ehc <- preds$next_log_ehc_z * orig_sd + orig_mean

Plot
class()
preds |>
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
  facet_wrap(Seafloor_bites_binary~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Step length"
  ) +
  theme_classic(base_size = 24)

full_data_strip$group_type <- as.factor(full_data_strip$group_type)
max(full_movement_data$step_length)
