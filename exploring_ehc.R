full_data %>%
  ggplot(aes(x = hit_rate, y = Seafloor_bites_binary)) +
  geom_point() +
  geom_smooth()+
  theme_minimal()
full_data %>%
  ggplot(aes(x = min_hit_distance, y = Seafloor_bites_binary)) +
  geom_point() +
  geom_smooth()

cor(full_data_strip$avg_norm_height, full_data_strip$avg_distance)

disttable <- full_data_strip %>%
  #filter(Seafloor_bites_binary == 1) %>%
  dplyr::select(Seafloor_bites_binary, ObsID, sec, min_hit_distance, log_ehc, log_ehc_height_only, hit_rate) %>%
  filter(min_hit_distance <0.22)

disttable %>%
  ggplot(aes(x = hit_rate, y = Seafloor_bites_binary)) +
  geom_point() +
  geom_smooth()

# Does height scale with distance within the same trial?
full_data_strip %>%
  group_by(Trial_ID) %>%
  summarise(
    cor_height_distance = cor(avg_norm_height, avg_distance, 
                              use = "complete.obs")
  ) %>%
  summary()

# Plot within trials
ggplot(full_data_strip,
       aes(avg_distance, avg_norm_height)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~ Trial_ID) +
  theme_classic() +
  labs(title = "Height vs distance within each trial")

full_data_strip %>%
  group_by(activity_state) %>%
  summarise(mean_hit_rate = mean(hit_rate, na.rm = TRUE))

full_data_strip$avg_norm_height_z <- scale(full_data_strip$avg_norm_height)
seafloor_bites_distance_model <- brm(
  formula = Seafloor_bites_binary ~
    avg_norm_height +
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

summary(seafloor_bites_distance_model)

test <- full_data_strip %>%
  filter(activity_state == "Foraging") %>%
  arrange(desc(avg_height)) %>%
  dplyr::select(Trial_ID, sec, frame_n, avg_height, avg_norm_height,
         min_hit_height, max_hit_height, avg_distance) %>%
  head(10)


plot(full_data_strip$hit_rate, full_data_strip$avg_distance)
plot(full_data_strip$hit_rate, full_data_strip$avg_norm_height)
plot(full_data_strip$avg_distance, full_data_strip$avg_norm_height)


corrplot::corrplot()

install.packages("GGally")
library(GGally)

cor_data <- full_data_strip %>%
  dplyr::select(
    hit_rate, 
    avg_distance, 
    avg_norm_height,
  ) %>%
  drop_na()


pairs(cor_data,
      pch = ".",
      col = rgb(0, 0, 0, 1),
      upper.panel = NULL)


library(lme4)

# Standardise predictors
full_data_strip$hit_rate_z <- as.numeric(scale(full_data_strip$hit_rate))
full_data_strip$avg_distance_z <- as.numeric(scale(full_data_strip$avg_distance))
full_data_strip$avg_norm_height_z <- as.numeric(scale(full_data_strip$avg_norm_height))

# Fit model
biting_model <- glmer(
  Seafloor_bites_binary ~
    avg_distance_z*subject*social_state +
    (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  data = full_data_strip,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl = list(maxfun = 2e5))
)

summary(biting_model)


library(lme4)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── Standardise predictors in both datasets ────────────────────────

full_data_strip$hit_rate_z        <- as.numeric(scale(full_data_strip$hit_rate))
full_data_strip$avg_distance_z    <- as.numeric(scale(full_data_strip$avg_distance))
full_data_strip$avg_norm_height_z <- as.numeric(scale(full_data_strip$avg_norm_height))

decision_data$hit_rate_z        <- as.numeric(scale(decision_data$hit_rate))
decision_data$avg_distance_z    <- as.numeric(scale(decision_data$avg_distance))
decision_data$avg_norm_height_z <- as.numeric(scale(decision_data$avg_norm_height))

# ── Fit biting models ──────────────────────────────────────────────

bite_hit <- glmer(
  Seafloor_bites_binary ~ hit_rate_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data = full_data_strip, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

bite_dist <- glmer(
  Seafloor_bites_binary ~ avg_distance_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data = full_data_strip, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

bite_height <- glmer(
  Seafloor_bites_binary ~ avg_norm_height_z * subject * social_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data = full_data_strip, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# ── Fit foraging initiation models ────────────────────────────────

init_hit <- glmer(
  next_state_different ~ hit_rate_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data = decision_data, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

init_dist <- glmer(
  next_state_different ~ avg_distance_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data = decision_data, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

init_height <- glmer(
  next_state_different ~ avg_norm_height_z * subject * next_social_state +
    time_in_state +
    (1 | Trial_ID) + (1 | site) + (1 | folder_date),
  data = decision_data, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# ── Plot function ──────────────────────────────────────────────────
plot_model_fit <- function(model, data, x_var, y_var,
                           x_label, y_label, social_var,
                           title = "") {
  
  # Prediction grid
  x_range <- seq(min(data[[x_var]], na.rm = TRUE),
                 max(data[[x_var]], na.rm = TRUE),
                 length.out = 100)
  
  social_levels <- sort(unique(data[[social_var]]))
  
  pred_grid <- expand.grid(
    x       = x_range,
    subject = c("Parrotfish", "Surgeonfish"),
    social  = social_levels,
    stringsAsFactors = FALSE
  )
  
  if ("time_in_state" %in% colnames(data)) {
    pred_grid$time_in_state <- median(data$time_in_state, na.rm = TRUE)
  }
  
  names(pred_grid)[names(pred_grid) == "x"]      <- x_var
  names(pred_grid)[names(pred_grid) == "social"] <- social_var
  
  pred_grid$Trial_ID    <- NA
  pred_grid$site        <- NA
  pred_grid$folder_date <- NA
  
  pred_grid$mnmu <- predict(model, newdata = pred_grid,
                            type = "response", re.form = NA)
  
  # Plot
  ggplot() +
    geom_jitter(
      data = data,
      aes(x     = .data[[x_var]],
          y     = .data[[y_var]],
          color = .data[[social_var]]),
      alpha  = 0.2,
      size   = 1,
      height = 0.01
    ) +
    geom_line(
      data = pred_grid,
      aes(x     = .data[[x_var]],
          y     = mnmu,
          color = .data[[social_var]]),
      linewidth = 1
    ) +
    facet_wrap(~ subject) +
    labs(
      x     = x_label,
      y     = y_label,
      color = "Social state",
      title = title
    ) +
    theme_classic(base_size = 24) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold"))
}

# ── Biting plots ───────────────────────────────────────────────────

p_bite_hit <- plot_model_fit(
  bite_hit, full_data_strip,
  x_var      = "hit_rate_z",
  y_var      = "Seafloor_bites_binary",
  x_label    = "Reef density",
  y_label    = "Probability of biting",
  social_var = "social_state",
  title      = "Reef density"
)

p_bite_dist <- plot_model_fit(
  bite_dist, full_data_strip,
  x_var      = "avg_distance_z",
  y_var      = "Seafloor_bites_binary",
  x_label    = "Distance to substrate",
  y_label    = "Probability of biting",
  social_var = "social_state",
  title      = "Distance to substrate"
)

p_bite_height <- plot_model_fit(
  bite_height, full_data_strip,
  x_var      = "avg_norm_height_z",
  y_var      = "Seafloor_bites_binary",
  x_label    = "Vertical heterogeneity",
  y_label    = "Probability of biting",
  social_var = "social_state",
  title      = "Vertical heterogeneity"
)

# ── Foraging initiation plots ──────────────────────────────────────

p_init_hit <- plot_model_fit(
  init_hit, decision_data,
  x_var      = "hit_rate_z",
  y_var      = "next_state_different",
  x_label    = "Reef density",
  y_label    = "Probability of initiating foraging",
  social_var = "next_social_state",
  title      = "Reef density"
)

p_init_dist <- plot_model_fit(
  init_dist, decision_data,
  x_var      = "avg_distance_z",
  y_var      = "next_state_different",
  x_label    = "Distance to substrate",
  y_label    = "Probability of initiating foraging",
  social_var = "next_social_state",
  title      = "Distance to substrate"
)

p_init_height <- plot_model_fit(
  init_height, decision_data,
  x_var      = "avg_norm_height_z",
  y_var      = "next_state_different",
  x_label    = "Vertical heterogeneity",
  y_label    = "Probability of initiating foraging",
  social_var = "next_social_state",
  title      = "Vertical heterogeneity"
)

# ── Combine ────────────────────────────────────────────────────────

(p_bite_dist | p_bite_height | p_bite_hit) /
  (p_init_dist | p_init_height | p_init_hit) +
  plot_annotation(tag_levels = "A")



full_data %>%
  ggplot()