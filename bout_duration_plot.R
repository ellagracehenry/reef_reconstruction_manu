Create the prediction grid

#bout duration
pred_newd <- foraging_bout_stats |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    bout_log_ehc_z = list(seq(min(bout_log_ehc_z), max(bout_log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(bout_log_ehc_z)

#bite rate
pred_newd <- foraging_bout_stats |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    bout_log_ehc_z = list(
      seq(
        min(bout_log_ehc_z),
        max(bout_log_ehc_z),
        length.out = 50
      )
    ),
    .groups = "drop"
  ) |>
  tidyr::unnest(bout_log_ehc_z) |>
  dplyr::mutate(
    bout_duration = 1
  )


Keep only subject × social_state combinations that actually occur in the data

observed_combinations <- foraging_bout_stats |>
  dplyr::distinct(subject, social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_state")
  )

Generate predictions

pmu <- posterior_epred(
  bite_count_mod,
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

orig_mean <- attr(foraging_bout_stats$bout_log_ehc_z, "scaled:center")
orig_sd <- attr(foraging_bout_stats$bout_log_ehc_z, "scaled:scale")

preds$bout_log_ehc <- preds$bout_log_ehc_z * orig_sd + orig_mean

Plot

preds |>
  ggplot() +
  geom_jitter(
    data = foraging_bout_stats,
    aes(
      x = bout_log_ehc,
      y = bout_bites,
      color = social_state
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  geom_ribbon(
    aes(
      x = bout_log_ehc,
      ymin = mulo95,
      ymax = muhi95,
      fill = social_state
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = bout_log_ehc,
      y = mnmu,
      color = social_state
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Total bites per bout"
  ) +
  theme_classic(base_size = 24)
