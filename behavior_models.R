
library(brms)
library(loo)
library(dplyr)
library(purrr)
library(stringr)

#BITES
full_data <- full_data 

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

full_data_strip <- full_data |>
  tidyr::drop_na(ehc_z, subject, modifier_1, modifier_2_z, Seafloor_bites_binary) 

plot(Seafloor_bites_binary ~ social_type, full_data_strip )

hist(full_data_strip$Seafloor_bites_binary)

hist(full_data$log_ehc_z)

levels(full_data_strip$social_state)
full_data_strip <- full_data_strip %>%
  mutate(
    social_type = case_when(
      social_state == "Social" & modifier_1 == "mixed family" ~ 2,
      social_state == "Social" & modifier_1 == "single family" ~ 1,
      social_state == "Asocial" ~ 0
    )
  )

full_data_strip$social_type <- as.factor(full_data_strip$social_type)

table(full_data_strip$social_type)
table(full_data_strip$social_state)
table(full_data_strip$modifier_1)

seafloor_bites_model_conspecifics <- brm(
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

pp_check(seafloor_bites_model_conspecifics, ndraws = 100)
summary(seafloor_bites_model_conspecifics)
posterior <- as.matrix(seafloor_bites_model_conspecifics)
mcmc_areas(posterior)
plot(seafloor_bites_model_conspecifics)
saveRDS(seafloor_bites_model_conspecifics, "binary_bites_mod.RData")


str(full_data_strip$folder_date)
full_data_strip <- full_data_strip %>%
  mutate(day_number = as.numeric(as.Date(as.character(folder_date))) - 
           as.numeric(as.Date("2023-06-08")))
full_data_strip$day_number


m_bayes <- brm(
  Seafloor_bites ~ log_ehc_z * subject * social_state + site + day_number + (1 | Trial_ID),
  data = full_data_strip,
  family = poisson(),
  prior = c(
    prior(normal(0, 0.5), class = b),           # fixed effects
    prior(normal(-0.5, 0.5), class = Intercept),   # intercept
    prior(exponential(2), class = sd)         # random effect SD
  ),
  #sample_prior = "only",
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 42
)

pp_check(m_bayes, type = "bars", ndraws = 100) + 
  xlim(-1, 10)
pp_check(m_bayes, type = "stat", stat = "mean", ndraws = 1000)
pp_check(m_bayes, type = "stat", stat = function(y) mean(y > 2), ndraws = 1000)
plot(m_bayes)
summary(m_bayes)


#Bite type
bite_type_data <- full_data_strip %>%
  filter(Water_column_bites + Seafloor_bites > 0)

bite_type_data <- bite_type_data |>
  dplyr::mutate(
    total_bites = Water_column_bites + Seafloor_bites,
    seafloor_prop = Seafloor_bites / total_bites
  )

bite_type_data$Seafloor_bites_binary
plot(seafloor_prop ~ log_ehc_z, data = bite_type_data)

bite_type_data$total_bites
bite_type_model <- brm(
  formula = Seafloor_bites_binary ~
    log_ehc_z*subject*social_type +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = bite_type_data,
  
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(normal(-2.5, 1.5), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  
  save_pars = save_pars(all = TRUE),
  
  sample_prior = "only",
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.95
  )
)

summary(bite_type_model)
pp_check(bite_type_model)

#FORAGING INITIATION

class(foraging_initiation_data$ehc_z)
class(foraging_initiation_data$foraging_initiation)
foraging_initiation_data$subject <- as.factor(foraging_initiation_data$subject)
foraging_initiation_data$social_state <- as.factor(foraging_initiation_data$social_state)
foraging_initiation_data$ehc_z <- scale(foraging_initiation_data$ehc)
foraging_initiation_data$log_ehc_z <- scale(log(foraging_initiation_data$ehc))
foraging_initiation_data$log_ehc <- log(foraging_initiation_data$ehc)

foraging_initiation_data <- foraging_initiation_data %>%
  mutate(subject = factor(
    subject,
    levels = c("Surgeonfish", "Parrotfish")
  )
  )

foraging_initiation_data <- foraging_initiation_data %>%
  mutate(subject = factor(
    subject,
    levels = c("Parrotfish","Surgeonfish")
  )
  )

levels(foraging_initiation_data$subject)


foraging_initiation_data_f <- foraging_initiation_data %>%
  filter(!ObsID %in% c(11,25))
foraging_initiation_data_f$ehc_z <- scale(foraging_initiation_data_f$ehc)


#Full
foraging_initiation_mod_full_log <- brm(formula = foraging_initiation ~
                        log_ehc_z*subject*social_state +
                        # ar(time = sec, gr = ObsID, p = 1, cov = TRUE) +
                        (1|ObsID) +
                        (1|site) +
                        (1|folder_date),
                      data = foraging_initiation_data,
                      chains = 4,
                      cores = 4,
                      iter = 4000,
                      warmup = 2000,
                      prior = c(
                        prior(normal(0, 1), class = "b"),
                        prior(normal(-2.5, 1.5), class = "Intercept"),
                        prior(exponential(1), class = "sd")
                      ),
                      save_pars = save_pars(all=TRUE),
                      #sample_prior = 'only',
                      family = bernoulli(link = "logit"),
                      control = list(adapt_delta = 0.98)
)
summary(foraging_initiation_mod_full_log)
pp_check(foraging_initiation_mod_full_log)

loo(foraging_initiation_mod_full, foraging_initiation_mod_full_log) #non log is 6.8 se diff worse

#Filtered 
foraging_initiation_mod_fil <- brm(formula = foraging_initiation ~
                                      log_ehc_z*subject*social_state +
                                      # ar(time = sec, gr = ObsID, p = 1, cov = TRUE) +
                                      (1|ObsID) +
                                      (1|site) +
                                      (1|folder_date),
                                    data = foraging_initiation_data_f,
                                    chains = 4,
                                    cores = 4,
                                    iter = 4000,
                                    warmup = 2000,
                                    prior = c(
                                      prior(normal(0, 1), class = "b"),
                                      prior(normal(-2.5, 1.5), class = "Intercept"),
                                      prior(exponential(1), class = "sd")
                                    ),
                                    save_pars = save_pars(all=TRUE),
                                    #sample_prior = 'only',
                                    family = bernoulli(link = "logit"),
                                    control = list(adapt_delta = 0.98)
)

summary(foraging_initiation_mod_fil)

#BOUT DURATION
bout_stats <- full_data %>%
  group_by(Trial_ID, bout, subject, folder_date, site, social_state
           #,modifiers
  ) %>%
  reframe(
    bout_type = first(activity_state),
    bout_duration = n(),                             # seconds in bout
    bout_bites = sum(Bites),                         # total bites taken
    bout_bite_rate = bout_bites / bout_duration,               # bites per second
    bout_ehc = mean(ehc, na.rm=TRUE), 
    bout_ehc_inverse_sd = mean(ehc_inverse_sd, na.rm=TRUE),
    bout_sd = mean(sd, na.rm=TRUE),
    bout_modifier_1 = first(modifier_1),
    bout_modifier_2 = mean(modifier_2, na.rm=TRUE),
    .groups = "drop"
  ) 

bout_stats <- bout_stats |>
  tidyr::drop_na(bout_bites, bout_duration, bout_ehc, social_state, bout_type)

foraging_bout_stats <- bout_stats %>%
  filter(bout_type == "Foraging") %>%
  filter(!Trial_ID %in% c(11,25))



foraging_bout_stats$bout_log_ehc_z <- scale(log(foraging_bout_stats$bout_ehc))
foraging_bout_stats$bout_log_ehc <- log(foraging_bout_stats$bout_ehc)

hist(bout_stats$bout_duration) #gamma or lognormal

foraging_bout_stats$social_state <- as.factor(foraging_bout_stats$social_state)
foraging_bout_stats$subject <- as.factor(foraging_bout_stats$subject)

bout_duration_mod_full <- brm(
  bout_duration ~ bout_log_ehc_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
    data = foraging_bout_stats,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    #prior(normal(0, 1), class = "b"),
    #prior(normal(0, 2), class = "Intercept"),
    #prior(exponential(1), class = "sd"),
    #prior(exponential(1), class = "shape")
    prior(normal(log(30), 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "sigma")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  #family = Gamma(link = "log"),
  family = lognormal(),
  control = list(adapt_delta = 0.98)
)

pp_check(bout_duration_mod, draw = 1000)
pp_check(bout_duration_mod, type = "hist")
pp_check(bout_duration_mod, type = "stat", stat = "mean")
plot(bout_duration)
pp_check(bout_duration_mod, type = "dens_overlay", n_draws = 2000) +
  ggplot2::scale_x_log10()

summary(bout_duration_mod)
pp_check(bout_duration_mod)

#BITE RATE
hist(bout_stats$bout_bites)

bite_rate_mod <- brm(
  bout_bites ~ bout_log_ehc_z*social_state*subject + offset(log(bout_duration)) +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = foraging_bout_stats,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
   # prior(normal(0, 1), class = "b"),
    #prior(normal(0, 2), class = "Intercept")
    prior(normal(log(3), 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "shape")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  family = negbinomial(),
  control = list(adapt_delta = 0.98)
)

summary(bite_rate_mod)
pp_check(bite_rate_mod)
var(foraging_bout_stats$bite_rate)
mean(foraging_bout_stats$bite_rate)
#if variance > mean then use negative binomial (better for overdispered)

#TOTAL BITES PER BOUT
bite_count_mod <- brm(
  bout_bites ~ bout_log_ehc_z*social_state*subject +
    (1|Trial_ID) +
    (1|site) +
    (1|folder_date),
  data = foraging_bout_stats,
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  prior = c(
    #prior(normal(0, 1), class = "b"), #poisson priors
    #prior(normal(0, 2), class = "Intercept")
    prior(normal(log(3), 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(exponential(1), class = "shape")
  ),
  save_pars = save_pars(all=TRUE),
  #sample_prior = 'only',
  family = negbinomial(),
  control = list(adapt_delta = 0.98)
)

summary(bite_count_mod)
pp_check(bite_count_mod)
var(foraging_bout_stats$bout_bites)
mean(foraging_bout_stats$bout_bites)
#if variance > mean then use negative binomial (better for overdispered)