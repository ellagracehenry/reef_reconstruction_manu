pred_newd <- full_data_strip |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

observed_combinations <- full_data_strip |>
  dplyr::distinct(subject, social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_state")
  )

pmu <- posterior_epred(
  activity_models$lESSs_lExS_lExSs_SxSs,
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

orig_mean <- attr(full_data_strip$log_ehc_z, "scaled:center")
orig_sd <- attr(full_data_strip$log_ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

Undo the log transformation

preds$ehc <- exp(preds$log_ehc)

Original-scale EHC for observed data

full_data_strip$ehc <- exp(full_data_strip$log_ehc)

Plot on the intuitive, original EHC scale

preds |>
  ggplot() +
  geom_jitter(
    data = full_data_strip,
    aes(
      x = ehc,
      y = bites_binary,
      color = social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = ehc,
      ymin = mulo95,
      ymax = muhi95,
      fill = social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = ehc,
      y = mnmu,
      color = social_state
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Probability of biting"
  ) +
  theme_classic(base_size = 24)
