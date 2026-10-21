#LIBRARIES
library(brms)
library(loo)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

#PREP
full_data <- full_data %>%
  filter(!Trial_ID %in% c(11,25))

full_data$ehc_z <- scale(full_data$ehc)
full_data$ehc_inverse_sd_z <- as.numeric(scale(full_data$ehc_inverse_sd))
full_data$sd_z <- as.numeric(scale(full_data$sd))
full_data$subject <- as.factor(full_data$subject)
full_data$modifier_1 <- as.factor(full_data$modifier_1)
full_data$modifier_2_z <- as.numeric(full_data$modifier_2)
full_data$Trial_ID <- as.factor(full_data$Trial_ID)
full_data$site <- as.factor(full_data$site)
full_data$folder_date <- as.factor(full_data$folder_date)
full_data$social_state <- as.factor(full_data$social_state)
full_data$ehc_squared_z <- as.numeric((full_data$ehc)^2)
full_data$log_ehc_z <- scale(log(full_data$ehc))
full_data$log_ehc <- log(full_data$ehc)
full_data$log_ehc_z_squared <- full_data$log_ehc_z^2

full_data <- full_data %>%
  mutate(
    social_type = case_when(
      social_state == "Social" & modifier_1 == "mixed family" ~ 2,
      social_state == "Social" & modifier_1 == "single family" ~ 1,
      social_state == "Asocial" ~ 0
    )
  )

full_data$social_type <- as.factor(full_data$social_type)

full_data_strip <- full_data |>
  tidyr::drop_na(ehc_z, subject, modifier_1, modifier_2_z, Seafloor_bites_binary) 

plot(Seafloor_bites_binary ~ social_type, full_data_strip )

hist(full_data_strip$Seafloor_bites_binary)

hist(full_data$log_ehc_z)

levels(full_data_strip$social_state)


full_data_strip <- full_data_strip %>%
  mutate(
      group_type = case_when(
      social_state == "Social" & modifier_2 == 2 ~ 1, #pair
      social_state == "Social" & modifier_2 > 2 ~ 3, #group
      social_state == "Asocial" ~ 0 #solo
    )
  )
table(full_data_strip$group_type)
full_data_strip$social_type <- as.factor(full_data_strip$social_type)



#MODELS
###1 - Seafloor bites ~ ehc, family, social vs asocial
class(full_data_strip$social_state)
seafloor_bites_model <- brm(
  formula = Seafloor_bites_binary ~
    log_ehc_z*subject*social_state +
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)
summary(seafloor_bites_model)
as_draws_df(seafloor_bites_model) |>
  mutate(surgeonfish_social = b_social_stateSocial + `b_subjectSurgeonfish:social_stateSocial`) |>
  summarise(
    mean = mean(surgeonfish_social),
    lower = quantile(surgeonfish_social, 0.025),
    upper = quantile(surgeonfish_social, 0.975)
  )
pp_check(seafloor_bites_model, ndraws=100)

saveRDS(seafloor_bites_model, file = "seafloor_bites_model.rds")

#fig
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
  seafloor_bites_model,
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
orig_mean <- attr(full_data_strip$log_ehc_z, "scaled:center")
orig_sd <- attr(full_data_strip$log_ehc_z, "scaled:scale")
preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean
bites_social_state <- preds |>
  ggplot() +
  geom_jitter(
    data = full_data_strip,
    aes(
      x = log_ehc,
      y = Seafloor_bites_binary,
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
    y = "Probability of biting"
  ) +
  theme_classic(base_size = 24)

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/bites_social_state.png", plot = bites_social_state, width = 10, height = 5)


#2 - Seafloor bites ~ ehc, family, solo vs conspecific group vs heterospecific group
table(full_data_strip$social_type)
class(full_data_strip$social_type)
bite_type_model <- brm(
  formula = Seafloor_bites_binary ~
    log_ehc_z*subject*social_type +
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)
saveRDS(bite_type_model, file = "seafloor_bites_social_type_model.rds")

summary(bite_type_model)
pp_check(bite_type_model, ndraws=100)

#fig
pred_newd <- full_data_strip |>
  dplyr::group_by(subject, social_type) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)
  
observed_combinations <- full_data_strip |>
  dplyr::distinct(subject, social_type)
pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_type")
  )
pmu <- posterior_epred(
  bite_type_model,
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
orig_mean <- attr(full_data_strip$log_ehc_z, "scaled:center")
orig_sd <- attr(full_data_strip$log_ehc_z, "scaled:scale")
preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

bites_social_type <- preds |>
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
  theme_classic(base_size = 24)

bites_social_type

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/bites_social_type.png", plot = bites_social_type, width = 10, height = 5)


#3 - Seafloor bites ~ ehc, family, solo vs pair vs group
full_data_strip$group_type <- as.factor(full_data_strip$group_type)
table(full_data_strip$group_type)

group_type_model <- brm(
  formula = Seafloor_bites_binary ~
    log_ehc_z*subject*group_type +
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)
summary(group_type_model)
saveRDS(group_type_model, file = "seafloor_bites_group_type_model.rds")

#fig
pred_newd <- full_data_strip |>
  dplyr::group_by(subject, group_type) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

observed_combinations <- full_data_strip |>
  dplyr::distinct(subject, group_type)
pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "group_type")
  )
pmu <- posterior_epred(
  group_type_model,
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
orig_mean <- attr(full_data_strip$log_ehc_z, "scaled:center")
orig_sd <- attr(full_data_strip$log_ehc_z, "scaled:scale")
preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

bites_group_type <- preds |>
  ggplot() +
  geom_jitter(
    data = full_data_strip,
    aes(
      x = log_ehc,
      y = Seafloor_bites_binary,
      color = group_type
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
      fill = group_type
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      x = log_ehc,
      y = mnmu,
      color = group_type
    ),
    linewidth = 1
  ) +
  facet_wrap(~subject) +
  labs(
    x = "Experienced Habitat Complexity",
    y = "Probability of biting"
  ) +
  theme_classic(base_size = 24)

bites_social_type

summary(group_type_model)
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/bites_group_type.png", plot = bites_group_type, width = 10, height = 5)



#Foraging initiation

#PREP FOR FORAGING INITIATION
bout_summary <- full_data %>%
  filter(activity_state == "Foraging") %>%
  group_by(ObsID, bout) %>%
  summarise(
    seafloor_bout = as.integer(
      any(Seafloor_bites > 0, na.rm = TRUE)
    ),
    water_column_bout = as.integer(
      any(Water_column_bites > 0, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Assign initiation to the preceding row and identify the upcoming foraging bout
foraging_initiation_data <- full_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  mutate(
    next_activity = lead(activity_state),
    next_bout = lead(bout),
    next_sec = lead(sec),
    
    foraging_initiation = as.integer(
      activity_state == "Traversing/swimming" &
        next_activity == "Foraging" &
        next_sec == sec + 1
    )
  ) %>%
  ungroup() %>%
  filter(
    activity_state  != "Traversing/swimming" |
      (activity_state == "Traversing/swimming" & !is.na(next_activity) & next_sec == sec + 1)
  ) %>%
  left_join(
    bout_summary %>%
      rename(
        next_bout = bout,
        next_seafloor_bout = seafloor_bout,
        next_water_column_bout = water_column_bout
      ),
    by = c("ObsID", "next_bout")
  ) %>%
  mutate(
    seafloor_bout_initiation = as.integer(
      foraging_initiation == 1 & !is.na(next_seafloor_bout) & next_seafloor_bout == 1
    ),
    water_column_bout_initiation = as.integer(
      foraging_initiation == 1 & !is.na(next_water_column_bout) & next_water_column_bout == 1
    )
  ) %>%
  select(-next_activity, -next_bout, -next_sec, -next_seafloor_bout,
         -next_water_column_bout)

traversing_foraging_initiation_data <- foraging_initiation_data %>%
  filter(activity_state == "Traversing/swimming")


traversing_foraging_initiation_data <- traversing_foraging_initiation_data |>
  tidyr::drop_na(log_ehc_z, subject, seafloor_bout_initiation, social_type) 

class(traversing_foraging_initiation_data$social_type)
table(traversing_foraging_initiation_data$social_type)
table(traversing_foraging_initiation_data$seafloor_bout_initiation)

#4 - Seafloor foraging initiation, social type
seafloor_initiations_model <- brm(
  formula = seafloor_bout_initiation ~
    log_ehc_z*subject*social_type +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = traversing_foraging_initiation_data,
  
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)

summary(seafloor_initiations_model)

class(traversing_foraging_initiation_data$social_state)


#5 - Seafloor foraging initiation, social state
seafloor_initiations_social_state_model <- brm(
  formula = seafloor_bout_initiation ~
    log_ehc_z*subject*social_state +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = traversing_foraging_initiation_data,
  
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)
summary(seafloor_initiations_social_state_model)

#DECISION-MAKING
#PREP FOR DECISION MAKING MODEL
#Column for if next state is the same, given the next state is 1 second ahead
#TRAVERSING -> FORAGING FUTURE
decision_data <- full_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  mutate(
    next_activity = lead(activity_state),
    next_bout = lead(bout),
    next_sec = lead(sec),
    next_log_ehc = lead(log_ehc),
    next_social_state = lead(social_state), 
    
    next_state_different = as.integer(
      activity_state != next_activity &
        next_sec == sec + 1
    )
  ) %>%
  ungroup() %>%
  filter(
    activity_state == "Traversing/swimming" & !is.na(next_activity) & next_sec == sec + 1 & !is.na(next_log_ehc)) 

decision_data$next_social_state <- as.factor(decision_data$next_social_state)
decision_data$next_log_ehc_z <- scale(decision_data$next_log_ehc)
decision_data$next_state_different <- as.numeric(decision_data$next_state_different)
decision_data$subject <- as.factor(decision_data$subject)

decision_data <- decision_data |>
  tidyr::drop_na(ehc_z, next_log_ehc_z, subject, next_social_state, social_state) 

decision_data %>%
    ggplot(aes(x = next_log_ehc_z, y = next_state_different)) +
    geom_point() +
    geom_smooth()

foraging_initiation_model <- brm(
  formula = next_state_different ~
    next_log_ehc_z*subject*next_social_state +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = decision_data,
  
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)

saveRDS(foraging_initiation_model, file = "foraging_initiation_model.rds")
summary(foraging_initiation_model)
pp_check(foraging_initiation_model)

#fig
pred_newd <- decision_data |>
  dplyr::group_by(subject, next_social_state) |>
  dplyr::summarise(
    next_log_ehc_z = list(seq(min(next_log_ehc_z), max(next_log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(next_log_ehc_z)

observed_combinations <- decision_data |>
  dplyr::distinct(subject, next_social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "next_social_state")
  )

pmu <- posterior_epred(
  foraging_initiation_model,
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

orig_mean <- attr(decision_data$next_log_ehc_z, "scaled:center")
orig_sd <- attr(decision_data$next_log_ehc_z, "scaled:scale")

preds$next_log_ehc <- preds$next_log_ehc_z * orig_sd + orig_mean

foraging_initiation_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = decision_data,
    aes(
      x = next_log_ehc,
      y = next_state_different,
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
    y = "Probability of initiating foraging"
  ) +
  theme_classic(base_size = 24)

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/foraging_initiation.png", plot = foraging_initiation_fig, width = 10, height = 5)

#TRAVERSING -> FORAGING FUTURE
decision_data <- full_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  mutate(
    next_activity = lead(activity_state),
    next_bout = lead(bout),
    next_sec = lead(sec),
    next_log_ehc = lead(log_ehc),
    next_social_state = lead(social_state), 
    
    next_state_different = as.integer(
      activity_state != next_activity &
        next_sec == sec + 1
    )
  ) %>%
  ungroup() %>%
  filter(
    activity_state == "Foraging" & !is.na(next_activity) & next_sec == sec + 1 & !is.na(next_log_ehc)) 

decision_data$next_social_state <- as.factor(decision_data$next_social_state)
decision_data$next_log_ehc_z <- scale(decision_data$next_log_ehc)
decision_data$next_state_different <- as.numeric(decision_data$next_state_different)
decision_data$subject <- as.factor(decision_data$subject)

decision_data %>%
  ggplot(aes(x = next_log_ehc_z, y = next_state_different)) +
  geom_point() +
  geom_smooth()

decision_data <- 

table(decision_data$next_state_different)

traversing_initiation_model <- brm(
  formula = next_state_different ~
    next_log_ehc_z*subject*next_social_state +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = decision_data,
  
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)

saveRDS(traversing_initiation_model, file = "traversing_initiation_model.rds")
summary(traversing_initiation_model)


#fig
pred_newd <- decision_data |>
  dplyr::group_by(subject, next_social_state) |>
  dplyr::summarise(
    next_log_ehc_z = list(seq(min(next_log_ehc_z), max(next_log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(next_log_ehc_z)

observed_combinations <- decision_data |>
  dplyr::distinct(subject, next_social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "next_social_state")
  )

pmu <- posterior_epred(
  traversing_initiation_model,
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

orig_mean <- attr(decision_data$next_log_ehc_z, "scaled:center")
orig_sd <- attr(decision_data$next_log_ehc_z, "scaled:scale")

preds$next_log_ehc <- preds$next_log_ehc_z * orig_sd + orig_mean

traversing_initiation_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = decision_data,
    aes(
      x = next_log_ehc,
      y = next_state_different,
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
    y = "Probability of initiating traversing"
  ) +
  theme_classic(base_size = 24)
traversing_initiation_fig
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/traversing_initiation.png", plot = traversing_initiation_fig, width = 10, height = 5)
summary(traversing_initiation_model)


#TRAVERSING -> FORAGING CURRENT
decision_data <- full_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  mutate(
    next_activity = lead(activity_state),
    next_bout = lead(bout),
    next_sec = lead(sec),
    next_log_ehc = lead(log_ehc),
    next_social_state = lead(social_state), 
    
    next_state_different = as.integer(
      activity_state != next_activity &
        next_sec == sec + 1
    )
  ) %>%
  ungroup() %>%
  filter(
    activity_state == "Traversing/swimming" & !is.na(next_activity) & next_sec == sec + 1 & !is.na(next_log_ehc)) 

decision_data$next_social_state <- as.factor(decision_data$next_social_state)
decision_data$next_log_ehc_z <- scale(decision_data$next_log_ehc)
decision_data$next_state_different <- as.numeric(decision_data$next_state_different)
decision_data$subject <- as.factor(decision_data$subject)

decision_data %>%
  ggplot(aes(x = log_ehc, y = next_state_different)) +
  geom_point() +
  geom_smooth()

foraging_initiation_current_state_model <- brm(
  formula = next_state_different ~
    log_ehc_z*subject*social_state +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = decision_data,
  
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)

saveRDS(foraging_initiation_current_state_model, file = "foraging_initiation_current_state_model.rds")
summary(foraging_initiation_current_state_model)
pp_check(foraging_initiation_current_state_model)

#fig
pred_newd <- decision_data |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

observed_combinations <- decision_data |>
  dplyr::distinct(subject, social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_state")
  )

pmu <- posterior_epred(
  foraging_initiation_current_state_model,
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

orig_mean <- attr(decision_data$log_ehc_z, "scaled:center")
orig_sd <- attr(decision_data$log_ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

foraging_initiation_current_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = decision_data,
    aes(
      x = log_ehc,
      y = next_state_different,
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
    y = "Probability of initiating foraging"
  ) +
  theme_classic(base_size = 24)
foraging_initiation_current_fig
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/foraging_initiation_current_fig.png", plot = foraging_initiation_current_fig, width = 10, height = 5)

#TRAVERSING -> FORAGING CURRENT
decision_data <- full_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  mutate(
    next_activity = lead(activity_state),
    next_bout = lead(bout),
    next_sec = lead(sec),
    next_log_ehc = lead(log_ehc),
    next_social_state = lead(social_state), 
    
    next_state_different = as.integer(
      activity_state != next_activity &
        next_sec == sec + 1
    )
  ) %>%
  ungroup() %>%
  filter(
    activity_state == "Foraging" & !is.na(next_activity) & next_sec == sec + 1 & !is.na(next_log_ehc)) 

decision_data$next_social_state <- as.factor(decision_data$next_social_state)
decision_data$next_log_ehc_z <- scale(decision_data$next_log_ehc)
decision_data$next_state_different <- as.numeric(decision_data$next_state_different)
decision_data$subject <- as.factor(decision_data$subject)

decision_data %>%
  ggplot(aes(x = next_log_ehc_z, y = next_state_different)) +
  geom_point() +
  geom_smooth()

table(decision_data$next_state_different)

traversing_initiation_current_state_model <- brm(
  formula = next_state_different ~
    log_ehc_z*subject*social_state +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = decision_data,
  
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
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)

summary(traversing_initiation_current_state_model)

#fig
pred_newd <- decision_data |>
  dplyr::group_by(subject, social_state) |>
  dplyr::summarise(
    log_ehc_z = list(seq(min(log_ehc_z), max(log_ehc_z), length.out = 50)),
    .groups = "drop"
  ) |>
  tidyr::unnest(log_ehc_z)

observed_combinations <- decision_data |>
  dplyr::distinct(subject, social_state)

pred_newd <- pred_newd |>
  dplyr::semi_join(
    observed_combinations,
    by = c("subject", "social_state")
  )

pmu <- posterior_epred(
  traversing_initiation_current_state_model,
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

orig_mean <- attr(decision_data$log_ehc_z, "scaled:center")
orig_sd <- attr(decision_data$log_ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

traversing_initiation_current_fig <- preds |>
  ggplot() +
  geom_jitter(
    data = decision_data,
    aes(
      x = log_ehc,
      y = next_state_different,
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
    y = "Probability of initiating traversing"
  ) +
  theme_classic(base_size = 24)
traversing_initiation_current_fig
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/traversing_initiation_current_fig.png", plot = traversing_initiation_current_fig, width = 10, height = 5)
