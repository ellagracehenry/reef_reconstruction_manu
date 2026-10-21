#Create the prediction grid

pred_newd <- full_data_strip |>
  dplyr::group_by(subject, social_type) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

#Keep only subject × social_state combinations that actually occur in the data

observed_combinations <- full_data_strip |>
  dplyr::distinct(subject, social_type)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_type")
  )

#Generate predictions

pmu <- posterior_epred(
  bite_type_model
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

#Undo standardisation

orig_mean <- attr(full_data_strip$log_ehc_z, "scaled:center")
orig_sd <- attr(full_data_strip$log_ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

#Plot

preds |>
  ggplot() +
  geom_jitter(
    data = full_data_strip,
    aes(
      x = log_ehc,
      y = Seafloor_bites_binary,
      color = social_type
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
      fill = social_type
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = log_ehc,
      y = mnmu,
      color = social_type
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Probability of biting"
  ) +
  theme_classic(base_size = 24) +
  labs

