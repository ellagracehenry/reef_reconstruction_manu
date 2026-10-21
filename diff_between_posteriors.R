ehc_ranges <- full_data_strip |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    min_ehc = min(log_ehc_z, na.rm = TRUE),
    max_ehc = max(log_ehc_z, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = social_state,
    values_from = c(min_ehc, max_ehc)
  ) |>
  dplyr::mutate(
    lower = pmax(min_ehc_Asocial, min_ehc_Social),
    upper = pmin(max_ehc_Asocial, max_ehc_Social)
  )


pred_newd <- ehc_ranges |>
  dplyr::select(subject, lower, upper) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    log_ehc_z = list(
      seq(lower, upper, length.out = 50)
    )
  ) |>
  dplyr::ungroup() |>
  tidyr::unnest(log_ehc_z) |>
  tidyr::expand_grid(
    social_state = factor(
      c("Asocial", "Social"),
      levels = levels(full_data_strip$social_state)
    )
  ) |>
  dplyr::select(subject, social_state, log_ehc_z)


#Posterior predicted probabilities

pmu <- posterior_epred(
  activity_mod_logehc_relvel_nonf,
  newdata = pred_newd,
  re_formula = NA
)

#Match Social and Asocial predictions at each EHC value

social_idx <- which(pred_newd$social_state == "Social")
asocial_idx <- which(pred_newd$social_state == "Asocial")

#P(biting | Social) - P(biting | Asocial)

p_diff <- pmu[, social_idx] - pmu[, asocial_idx]

#Summarise posterior contrast

contrast_preds <- pred_newd[social_idx, ] |>
  dplyr::select(log_ehc_z, subject) |>
  dplyr::mutate(
    mean_diff = apply(p_diff, 2, mean),
    lo95 = apply(p_diff, 2, quantile, probs = 0.025),
    hi95 = apply(p_diff, 2, quantile, probs = 0.975)
  )

#Undo standardisation

orig_mean <- attr(full_data_strip$log_ehc_z, "scaled:center")
orig_sd <- attr(full_data_strip$log_ehc_z, "scaled:scale")

contrast_preds$log_ehc <-
  contrast_preds$log_ehc_z * orig_sd + orig_mean

#Undo log transformation

#contrast_preds$ehc <- exp(contrast_preds$log_ehc)

#Plot on original/unlogged EHC scale

contrast_preds |>
  ggplot(aes(x = log_ehc, y = mean_diff)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_ribbon(
    aes(ymin = lo95, ymax = hi95),
    alpha = 0.2
  ) +
  geom_line(
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Difference in probability of biting\nSocial − Asocial"
  ) +
  theme_classic(base_size = 24)

