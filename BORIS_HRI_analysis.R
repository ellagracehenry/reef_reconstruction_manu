#Packages
library(rstanarm)
library(brms)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(DHARMa)
library(tidyr)
library(loo)

#Functions
read_and_clean <- function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c")) %>%
    janitor::clean_names() %>%
    dplyr::mutate(
      source_file = basename(f),
      n = dplyr::row_number()
    )
}

hpdi <- function (samp, prob = 0.95) {
  vals <- sort(samp)
  nsamp <- length(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- which.min(vals[init + gap,drop=FALSE] - vals[init, drop=FALSE])
  ans <- cbind(lower=vals[inds], upper=vals[inds + gap])
  return(ans)
}

#read in boris df
boris23_timeline <- read.csv("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/BORIS/boris23_timeline_group_size.csv")
boris23_timeline$ObsName <- paste(boris23_timeline$folder, boris23_timeline$Trial_ID, sep = "_")
boris23_timeline$ObsID <- boris23_timeline$Trial_ID

#read i n conversions
conversions <- read.csv("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/BORIS/2023FF_Assignments.csv")
conversions <- conversions[ , c(3, 26, 27, 28)]
conversions$ObsName <- conversions$ObsID
conversions$ObsID <- NULL

#RAYCASTING output - read HRI
setwd("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/ehc")
csv_files <- list.files(pattern = "\\.csv$")

hri_data <- lapply(csv_files, read_and_clean) %>% bind_rows()
hri_data$ObsID <- as.numeric(str_extract_all(str_sub(hri_data$source_file, start = 11, end = 13), "\\d+"))

length(unique(hri_data$ObsID))

#MVT output - read in frame numbers and reconstruction/reprojection error
setwd("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/BORIS/frame_idx")
csv_files1 <- list.files(pattern = "\\.csv$")
frame_idx_data <- lapply(csv_files1, read_and_clean) %>% bind_rows()
frame_idx_data$ObsID <- as.numeric(str_extract_all(str_sub(frame_idx_data$source_file, start = 11, end = 13), "\\d+"))
frame_idx_data$frame_n <- as.numeric(frame_idx_data$frame)

unique(frame_idx_data$ObsID)

#Correct frame number for BORIS
boris23_timeline <- left_join(boris23_timeline, conversions, by = "ObsName")
boris23_timeline$frame_n <- (ceiling((((boris23_timeline$sec*3) + boris23_timeline$edit) - boris23_timeline$start_image_shift)/3)*3) - 2
#Foraging bouts manipulation
boris23_timeline <- boris23_timeline %>%
  arrange(Trial_ID, sec) %>%
  group_by(Trial_ID, subject) %>%
  mutate(
    bout = cumsum(activity_state != lag(activity_state, default = first(activity_state))))%>%
  ungroup()

#Correct frame number for HRI
hri_data <- left_join(frame_idx_data, hri_data, by = c("ObsID", "n"))

hri_data$reconstruction_error <- as.numeric(hri_data$reconstruction_error)
hri_data$reprojection_error <- as.numeric(hri_data$reprojection_error)
hri_data$min_hit_distance <- as.numeric(hri_data$min_hit_distance)
hri_data$max_hit_distance <- as.numeric(hri_data$max_hit_distance)
hri_data$min_hit_height <- as.numeric(hri_data$min_hit_height)
hri_data$max_hit_height <- as.numeric(hri_data$max_hit_height)
hri_data$avg_distance <- as.numeric(hri_data$avg_distance)
hri_data$avg_norm_height <- as.numeric(hri_data$avg_norm_height)

hri_data_f <- hri_data %>%
  filter(reprojection_error < 10) %>% #DECIDE ON RERPROJEC ERR CUT OFF 0.0246423
  filter(min_hit_distance > 0.01) #DECIDE ON DIST CUT OFF. 0.0464548

hri_data_f$ehc_constant_cv <- as.numeric(hri_data_f$ehc_constant_cv)
hri_data_f$ehc_constant_sd <- as.numeric(hri_data_f$ehc_constant_sd)
hri_data_f$ehc_inverse_cv <- as.numeric(hri_data_f$ehc_inverse_cv)
hri_data_f$ehc_inverse_sd <- as.numeric(hri_data_f$ehc_inverse_sd)
hri_data_f$ehc <- as.numeric(hri_data_f$ehc)
hri_data_f$ehc_add_sd <- as.numeric(hri_data_f$ehc_add_sd)
hri_data_f$ehc_sd <- as.numeric(hri_data_f$ehc_sd)
hri_data_f$cv <- as.numeric(hri_data_f$cv)
hri_data_f$sd <- as.numeric(hri_data_f$sd)




full_data <- left_join(boris23_timeline, hri_data_f, by = c("ObsID","frame_n"))
stripped_data <- full_data[, c("ObsID", "ObsName", "ehc", "ehc_add_sd", "ehc_sd","ehc_constant_cv","ehc_constant_sd","ehc_inverse_cv","ehc_inverse_sd","cv", "sd","sec","frame_n", "frame.y", "activity_state", "habitat_state", "social_state", "observation_date", "subject", "site", "folder_date", "folder_initials", "Bites","Sand.rubbing", "Departure", "Aggressive.on.focal", "Other.behavior", "Change.in.Focal.Fish", "Bite_classification", "modifier_1", "modifier_2", "hit_rate", "min_hit_distance", "max_hit_distance", "avg_distance", "avg_height", "avg_norm_height", "avg_distance_weighted_heights")]
stripped_data <- stripped_data %>% filter(activity_state %in% c("Traversing/swimming", "Foraging"))
stripped_data$activity_state <- as.factor(stripped_data$activity_state)

stripped_data <- stripped_data %>% mutate(modifier_1 = replace_na(modifier_1, "single family"))
stripped_data <- stripped_data %>% mutate(modifier_2 = replace_na(modifier_2, 0))
stripped_data$modifier_2 <- stripped_data$modifier_2 + 1


stripped_data %>%
  ggplot(aes(x = avg_distance, y = ehc_constant_cv)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = avg_norm_height, y = ehc_constant_cv)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = min_hit_distance, y = ehc_constant_cv)) +
  geom_point() +
  theme_classic()

stripped_data %>%
  ggplot(aes(x = avg_distance, y = ehc_constant_sd)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = avg_norm_height, y = ehc_constant_sd)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = min_hit_distance, y = ehc_constant_sd)) +
  geom_point() +
  theme_classic()

hri_data %>%
  ggplot(aes(x = as.numeric(avg_distance), y = as.numeric(ehc_inverse_cv))) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = avg_norm_height, y = ehc_inverse_cv)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = min_hit_distance, y = ehc_inverse_cv)) +
  geom_point() +
  theme_classic()


stripped_data %>%
  ggplot(aes(x = avg_distance, y = ehc_inverse_sd)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = avg_norm_height, y = ehc_inverse_sd)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = min_hit_distance, y = ehc_inverse_sd)) +
  geom_point() +
  theme_classic()

stripped_data %>%
  ggplot(aes(x = avg_distance, y = cv)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = avg_norm_height, y = cv)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = min_hit_distance, y = cv)) +
  geom_point() +
  theme_classic()

stripped_data %>%
  ggplot(aes(x = avg_distance, y = sd)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = avg_norm_height, y = sd)) +
  geom_point() +
  theme_classic()
stripped_data %>%
  ggplot(aes(x = min_hit_distance, y = sd)) +
  geom_point() +
  theme_classic()


stripped_data %>%
  ggplot(aes(x = as.numeric(cv), y = as.numeric(ehc))) +
  geom_point() +
  theme_classic()


plot_data <- stripped_data %>%
  select(
    avg_distance,
    avg_norm_height,
    min_hit_distance,
    ehc,
    ehc_add_sd,
    ehc_sd,
    ehc_constant_cv,
    ehc_constant_sd,
    ehc_inverse_cv,
    ehc_inverse_sd,
    cv,
    sd
  ) %>%
  pivot_longer(
    cols = c(
      ehc,
      ehc_add_sd,
      ehc_sd,
      ehc_constant_cv,
      ehc_constant_sd,
      ehc_inverse_cv,
      ehc_inverse_sd,
      cv,
      sd
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  pivot_longer(
    cols = c(
      avg_distance,
      avg_norm_height,
      min_hit_distance
    ),
    names_to = "predictor",
    values_to = "x"
  )

ggplot(plot_data, aes(x = x, y = as.numeric(value))) +
  geom_point(alpha = 0.5) +
  facet_grid(
    metric ~ predictor,
    scales = "free"
  ) +
  theme_classic(base_size = 9.5) +
  labs(
    x = NULL,
    y = NULL
  ) 

ggplot(plot_data, aes(x = as.numeric(value))) +
  geom_histogram() +
  facet_wrap(
  ~metric, scales="free"
  ) +
  theme_classic(base_size = 9.5) +
  labs(
    x = "value"
  ) 

#Distribution figure
#EHC value across trials
stripped_data %>%
  ggplot(aes(x = ehc_inverse_sd, fill = ObsID)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  labs(
    x = "Experienced Habitat Complexity"
  ) +
  theme_classic()


#ehc_constant_cv
#ehc_constant_sd
#ehc_inverse_cv
#ehc_inverse_sd
#min_hit_height
#max_hit_height
#min_hit_distance
#CV
#SD
#EHC
#EHC_add_sd
#EHC_sd

#hri_data <-  subset(hri_data, !(ObsID == 29 & hri_data$frame.y %in% c(142, 25, 212)))
#hri_data <-  subset(hri_data, !(ObsID == 40 & hri_data$frame.y %in% c(39)))
#hri_data <-  subset(hri_data, !(ObsID == 9 & hri_data$frame.y %in% c(78,77,75,97,90,40,138,165,164,162,160,159,229, 307)))
#hri_data <-  subset(hri_data, !(ObsID == 23 & hri_data$frame.y %in% c(1,261,263,535,23,190,187,188,189,190,183,273,357,360)))
#hri_data <-  subset(hri_data, !(ObsID == 61 & hri_data$frame.y %in% c(36, 45)))
#hri_data <-  subset(hri_data, !(ObsID == 41 & hri_data$frame.y %in% c(231,251,443,394,390)))
#hri_data <-  subset(hri_data, !(ObsID == 10 & hri_data$frame.y %in% c(63)))
#hri_data <-  subset(hri_data, !(ObsID == 66 & hri_data$frame.y %in% c(196)))
#hri_data <-  subset(hri_data, !(ObsID == 11 & hri_data$frame.y %in% c(115)))
#hri_data <-  subset(hri_data, !(ObsID == 39 & hri_data$frame.y %in% c(102, 103, 159, 215, 100)))
#hri_data <-  subset(hri_data, !(ObsID == 7 & hri_data$frame.y %in% c(136)))
#hri_data <-  subset(hri_data, !(ObsID == 27 & hri_data$frame.y %in% c(12, 153, 146, 90, 73, 62)))
#hri_data <-  subset(hri_data, !(ObsID == 64 & hri_data$frame.y %in% c(0)))
#hri_data <-  subset(hri_data, !(ObsID == 65 & hri_data$frame.y %in% c(252,40,11,112)))
#high HRI points are legit - EH30, 195, EH40, 162

#Bout dataframe
bout_stats <- full_data %>%
  group_by(Trial_ID, bout, subject, folder_date, site, social_state
           #,modifiers
  ) %>%
  summarize(
    bout_type = first(activity_state),
    bout_duration = n(),                             # seconds in bout
    bout_bites = sum(Bites),                         # total bites taken
    bout_bite_rate = bout_bites / bout_duration,               # bites per second
    bout_hri = mean(ehc_inverse_sd), 
    .groups = "drop"
  )

#Now scale the predictor variables for stability
bout_stats$bout_ehc_inverse_sd_z <- scale(bout_stats$bout_hri)
stripped_data$ehc_inverse_sd_z <- scale(stripped_data$ehc_inverse_sd)
stripped_data$activity_state <- as.factor(stripped_data$activity_state)
stripped_data$social_state <- as.factor(stripped_data$social_state)

#Model 1 - Foraging initiation. Now manipulate the data to indicate when a foraging initiation event occured
foraging_initiation_data <- stripped_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  mutate(
    activity_state_numeric = ifelse(
      activity_state == "Foraging", 0, 1
    ),
    activity_state_numeric_diff =
      activity_state_numeric - lag(activity_state_numeric),
    foraging_initiation =
      ifelse(activity_state_numeric_diff == -1, 1, 0)
  ) %>%
  filter(!(activity_state == "Foraging" & foraging_initiation == 0)) %>% #only keep traversing steps and the foraging step they started foraging on
  ungroup()

foraging_initiation_data$ehc_inverse_sd_z <- scale(foraging_initiation_data$ehc_inverse_sd)
hist(foraging_initiation_data$ehc_inverse_sd)
hist(foraging_initiation_data$ehc_inverse_sd_z)
foraging_initiation_data$modifier_2_z <- scale(foraging_initiation_data$modifier_2)
foraging_initiation_data$ehc_z <- scale(foraging_initiation_data$ehc)

foraging_initiation_data <- foraging_initiation_data |>
  tidyr::drop_na(foraging_initiation, ehc_inverse_sd_z, social_state, modifier_1)

table(foraging_initiation_data$modifier_2)

foraging_initiation_data <- foraging_initiation_data %>% mutate(
  group_size_binned = case_when(
    modifier_2 == 1 ~ "Solo",
    modifier_2 %in% c(2,3) ~ "Small",
    modifier_2 %in% c(4,5) ~ "Medium",
    modifier_2 %in% c(6,7) ~ "Large"
  ))

foraging_initiation_data$group_size_binned <- factor(group_size_binned, levels = c("Solo", "Small", "Medium", "Large"))

stripped_data <- stripped_data %>%
  mutate(
    group_size_binned = factor(
      group_size_binned,
      levels = c("Solo", "Small", "Medium", "Large")
    )
  )

foraging_initiation_data <- foraging_initiation_data %>%
  mutate(
    group_size_binned = factor(
      group_size_binned,
      levels = c("Solo", "Small", "Medium", "Large")
    ),
    
    modifier_1 = factor(
      modifier_1,
      levels = c("mixed family", "single family")
    ),
    
    subject = factor(
      subject,
      levels = c("Parrotfish", "Surgeonfish")
    )
  )

foraging_initiation_data$sd_z <- scale(foraging_initiation_data$sd)
hist(foraging_initiation_data$sd)
#Prior predictive checks on full model
activity_mod_6 <- brm(formula = foraging_initiation ~
                      sd_z +
                      group_size_binned +
                      modifier_1 + 
                     # subject +
                    #  ehc_z:group_size_binned +
                    #  ehc_z:modifier_1 +
                      sd_z:subject + 
                      sd_z:group_size_binned +
                   #   group_size_binned:modifier_1 +
                    #  group_size_binned:subject +
                    #  modifier_1:subject +
                      #ehc_z:group_size_binned:modifier_1 +
                   #   ehc_z:group_size_binned:subject +
                   #   group_size_binned:subject:modifier_1 +
                    #  modifier_1:subject:ehc_z +
                    #  ehc_z:modifier_1:subject +
                      #ehc_z:subject:group_size_binned:modifier_1 +
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
                      prior(normal(-2.5, 1.5), class = "Intercept")
                      #prior(exponential(1), class = "sd")
                    ),
                    save_pars = save_pars(all=TRUE),
                    #sample_prior = 'only',
                    family = bernoulli(link = "logit"),
                    control = list(adapt_delta = 0.98)
)

summary(activity_mod_6)
loo(activity_mod, activity_mod_2, activity_mod_3, activity_mod_4, activity_mod_5, activity_mod_6)

#PRIOR CHECKS
prior_pred <- posterior_predict(activity_mod)
dim(prior_pred)
prior_rate <- apply(prior_pred, 1, mean)

quantile(
  prior_rate,
  probs = c(.01, .05, .25, .5, .75, .95, .99)
)
#real initiation rate in range

hist(
  prior_rate,
  breaks = 50,
  main = "Prior predictive distribution",
  xlab = "Proportion initiating"
)

summary(activity_mod)
plot(activity_mod)

pp_check(activity_mod, type = "bars")
pp_check(activity_mod, type="dens_overlay",ndraws=500)

#Nested models
#Single factor models
activity_mod_E <- brm(formula = foraging_initiation ~
                        ehc_z +
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
                        prior(normal(-2.5, 1.5), class = "Intercept")
                        #prior(exponential(1), class = "sd")
                      ),
                      save_pars = save_pars(all=TRUE),
                      #sample_prior = 'only',
                      family = bernoulli(link = "logit"),
                      control = list(adapt_delta = 0.98)
)

activity_mod_G <- brm(formula = foraging_initiation ~
                           group_size_binned +
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
                           prior(normal(-2.5, 1.5), class = "Intercept")
                           #prior(exponential(1), class = "sd")
                         ),
                         save_pars = save_pars(all=TRUE),
                         #sample_prior = 'only',
                         family = bernoulli(link = "logit"),
                         control = list(adapt_delta = 0.98)
)

activity_mod_M <- brm(formula = foraging_initiation ~
                           modifier_1 + 
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
                           prior(normal(-2.5, 1.5), class = "Intercept")
                           #prior(exponential(1), class = "sd")
                         ),
                         save_pars = save_pars(all=TRUE),
                         #sample_prior = 'only',
                         family = bernoulli(link = "logit"),
                         control = list(adapt_delta = 0.98)
)

activity_mod_S <- brm(formula = foraging_initiation ~
                           subject +
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
                           prior(normal(-2.5, 1.5), class = "Intercept")
                           #prior(exponential(1), class = "sd")
                         ),
                         save_pars = save_pars(all=TRUE),
                         #sample_prior = 'only',
                         family = bernoulli(link = "logit"),
                         control = list(adapt_delta = 0.98)
)

#Two factor models
activity_mod_EG <- brm(formula = foraging_initiation ~
                           ehc_z +
                           group_size_binned +
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
                           prior(normal(-2.5, 1.5), class = "Intercept")
                           #prior(exponential(1), class = "sd")
                         ),
                         save_pars = save_pars(all=TRUE),
                         #sample_prior = 'only',
                         family = bernoulli(link = "logit"),
                         control = list(adapt_delta = 0.98)
)

activity_mod_GM <- brm(formula = foraging_initiation ~
                           group_size_binned +
                           modifier_1 + 
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
                           prior(normal(-2.5, 1.5), class = "Intercept")
                           #prior(exponential(1), class = "sd")
                         ),
                         save_pars = save_pars(all=TRUE),
                         #sample_prior = 'only',
                         family = bernoulli(link = "logit"),
                         control = list(adapt_delta = 0.98)
)

activity_mod_MS <- brm(formula = foraging_initiation ~
                         modifier_1 +
                         subject +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_SE <- brm(formula = foraging_initiation ~
                         ehc_z +
                         subject +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_SG <- brm(formula = foraging_initiation ~
                         group_size_binned +
                         subject +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_EM <- brm(formula = foraging_initiation ~
                         ehc_z +
                         modifier_1 +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

#Two factor models with interactions
activity_mod_ExG <- brm(formula = foraging_initiation ~
                         ehc_z +
                         group_size_binned +
                         ehc_z:group_size_binned +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_GxM <- brm(formula = foraging_initiation ~
                         group_size_binned +
                         modifier_1 + 
                         group_size_binned:modifier_1 +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_MxS <- brm(formula = foraging_initiation ~
                         modifier_1 +
                         subject +
                         modifier_1:subject + 
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_SxE <- brm(formula = foraging_initiation ~
                         ehc_z +
                         subject +
                         ehc_z:subject +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_SxG <- brm(formula = foraging_initiation ~
                         group_size_binned +
                         subject +
                         group_size_binned:subject +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_ExM <- brm(formula = foraging_initiation ~
                         ehc_z +
                         modifier_1 +
                         ehc_z:modifier_1 +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

################################

#Three factor models
activity_mod_EMG <- brm(formula = foraging_initiation ~
                         ehc_z +
                         modifier_1 +
                         group_size_binned +
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
                         prior(normal(-2.5, 1.5), class = "Intercept")
                         #prior(exponential(1), class = "sd")
                       ),
                       save_pars = save_pars(all=TRUE),
                       #sample_prior = 'only',
                       family = bernoulli(link = "logit"),
                       control = list(adapt_delta = 0.98)
)

activity_mod_EMS <- brm(formula = foraging_initiation ~
                          ehc_z +
                          modifier_1 +
                          group_size_binned +
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

activity_mod_EGS <- brm(formula = foraging_initiation ~
                          ehc_z +
                          group_size_binned +
                          subject +
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

activity_mod_MSG <- brm(formula = foraging_initiation ~
                          group_size_binned +
                          subject +
                          modifier_1 + 
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

#Three factor models with single two-way interactions
activity_mod_EMG_ExM <- brm(formula = foraging_initiation ~
                          ehc_z +
                          modifier_1 +
                          group_size_binned +
                          ehc_z:modifier_1 +
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

activity_mod_EMG_GxM <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              group_size_binned:modifier_1 +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EMG_ExG <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              ehc_z:group_size_binned +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)


activity_mod_EMS_ExM <- brm(formula = foraging_initiation ~
                          ehc_z +
                          modifier_1 +
                          group_size_binned +
                          ehc_z:modifier_1 +
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

activity_mod_EMS_ExG <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              ehc_z:group_size_binned +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EMS_ExM <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              group_size_binned:modifier_1 +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)


activity_mod_EGS_ExG <- brm(formula = foraging_initiation ~
                          ehc_z +
                          group_size_binned +
                          subject +
                          ehc_z:group_size_binned + 
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

activity_mod_EGS_ExS <- brm(formula = foraging_initiation ~
                          ehc_z +
                          group_size_binned +
                          subject +
                          ehc_z:subject + 
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

activity_mod_EGS_GxS <- brm(formula = foraging_initiation ~
                              ehc_z +
                              group_size_binned +
                              subject +
                              group_size_binned:subject + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_MSG_GxS <- brm(formula = foraging_initiation ~
                          group_size_binned +
                          subject +
                          modifier_1 + 
                          group_size_binned:subject + 
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
                          prior(normal(-2.5, 1.5), class = "Intercept")
                          #prior(exponential(1), class = "sd")
                        ),
                        save_pars = save_pars(all=TRUE),
                        #sample_prior = 'only',
                        family = bernoulli(link = "logit"),
                        control = list(adapt_delta = 0.98)
)

activity_mod_MSG_MxS <- brm(formula = foraging_initiation ~
                              group_size_binned +
                              subject +
                              modifier_1 + 
                              modifier_1:subject + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_MSG_GxM <- brm(formula = foraging_initiation ~
                              group_size_binned +
                              subject +
                              modifier_1 + 
                              group_size_binned:modifier_1 + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

###########

#Three factor model with multiple two way interactions
activity_mod_EMG_ExM <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              ehc_z:modifier_1 +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EMG_ExM <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              ehc_z:modifier_1 +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EMG_GxM <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              group_size_binned:modifier_1 +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EMG_ExG <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              ehc_z:group_size_binned +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)


activity_mod_EMS_ExM <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              ehc_z:modifier_1 +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EMS_ExG <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              ehc_z:group_size_binned +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EMS_ExM <- brm(formula = foraging_initiation ~
                              ehc_z +
                              modifier_1 +
                              group_size_binned +
                              group_size_binned:modifier_1 +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)


activity_mod_EGS_ExG <- brm(formula = foraging_initiation ~
                              ehc_z +
                              group_size_binned +
                              subject +
                              ehc_z:group_size_binned + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EGS_ExS <- brm(formula = foraging_initiation ~
                              ehc_z +
                              group_size_binned +
                              subject +
                              ehc_z:subject + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_EGS_GxS <- brm(formula = foraging_initiation ~
                              ehc_z +
                              group_size_binned +
                              subject +
                              group_size_binned:subject + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_MSG_GxS <- brm(formula = foraging_initiation ~
                              group_size_binned +
                              subject +
                              modifier_1 + 
                              group_size_binned:subject + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_MSG_MxS <- brm(formula = foraging_initiation ~
                              group_size_binned +
                              subject +
                              modifier_1 + 
                              modifier_1:subject + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_MSG_GxM <- brm(formula = foraging_initiation ~
                              group_size_binned +
                              subject +
                              modifier_1 + 
                              group_size_binned:modifier_1 + 
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_ESsSExSsExSSxSExSxSs <- brm(formula = foraging_initiation ~
                      ehc_inverse_sd_z +
                      social_state +
                      subject +
                      ehc_inverse_sd_z:social_state +
                      ehc_inverse_sd_z:subject + 
                      social_state:subject +
                      ehc_inverse_sd_z:subject:social_state +
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
                      prior(normal(-2.5, 1.5), class = "Intercept")
                      #prior(exponential(1), class = "sd")
                    ),
                    save_pars = save_pars(all=TRUE),
                    #sample_prior = 'only',
                    family = bernoulli(link = "logit"),
                    control = list(adapt_delta = 0.98)
)

activity_mod_ESsSExSsExSSxS <- brm(formula = foraging_initiation ~
                      ehc_inverse_sd_z +
                      social_state +
                      subject +
                      ehc_inverse_sd_z:social_state +
                      ehc_inverse_sd_z:subject + 
                      social_state:subject +
                      #ehc_inverse_sd_z:subject:social_state +
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
                      prior(normal(-2.5, 1.5), class = "Intercept")
                      #prior(exponential(1), class = "sd")
                    ),
                    save_pars = save_pars(all=TRUE),
                    #sample_prior = 'only',
                    family = bernoulli(link = "logit"),
                    control = list(adapt_delta = 0.98)
)
activity_mod_ESsSExSsExS <- brm(formula = foraging_initiation ~
                                     ehc_inverse_sd_z +
                                     social_state +
                                     subject +
                                     ehc_inverse_sd_z:social_state +
                                     ehc_inverse_sd_z:subject + 
                                     #social_state:subject +
                                     #ehc_inverse_sd_z:subject:social_state +
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
                                     prior(normal(-2.5, 1.5), class = "Intercept")
                                     #prior(exponential(1), class = "sd")
                                   ),
                                   save_pars = save_pars(all=TRUE),
                                   #sample_prior = 'only',
                                   family = bernoulli(link = "logit"),
                                   control = list(adapt_delta = 0.98)
)

activity_mod_ESsSExSsSsxS <- brm(formula = foraging_initiation ~
                                  ehc_inverse_sd_z +
                                  social_state +
                                  subject +
                                  ehc_inverse_sd_z:social_state +
                                  #ehc_inverse_sd_z:subject + 
                                  social_state:subject +
                                  #ehc_inverse_sd_z:subject:social_state +
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
                                  prior(normal(-2.5, 1.5), class = "Intercept")
                                  #prior(exponential(1), class = "sd")
                                ),
                                save_pars = save_pars(all=TRUE),
                                #sample_prior = 'only',
                                family = bernoulli(link = "logit"),
                                control = list(adapt_delta = 0.98)
)

activity_mod_ESsSExSSsxS <- brm(formula = foraging_initiation ~
                                  ehc_inverse_sd_z +
                                  social_state +
                                  subject +
                                  #ehc_inverse_sd_z:social_state +
                                  ehc_inverse_sd_z:subject + 
                                  social_state:subject +
                                  #ehc_inverse_sd_z:subject:social_state +
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
                                  prior(normal(-2.5, 1.5), class = "Intercept")
                                  #prior(exponential(1), class = "sd")
                                ),
                                save_pars = save_pars(all=TRUE),
                                #sample_prior = 'only',
                                family = bernoulli(link = "logit"),
                                control = list(adapt_delta = 0.98)
)

activity_mod_ESsSExSs <- brm(formula = foraging_initiation ~
                              ehc_inverse_sd_z +
                              social_state +
                              subject +
                              ehc_inverse_sd_z:social_state +
                              #ehc_inverse_sd_z:subject + 
                              #social_state:subject +
                              #ehc_inverse_sd_z:subject:social_state +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)


activity_mod_ESsSExS <- brm(formula = foraging_initiation ~
                              ehc_inverse_sd_z +
                              social_state +
                              subject +
                              #ehc_inverse_sd_z:social_state +
                              ehc_inverse_sd_z:subject + 
                              #social_state:subject +
                              #ehc_inverse_sd_z:subject:social_state +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)


activity_mod_ESsSSsxS <- brm(formula = foraging_initiation ~
                              ehc_inverse_sd_z +
                              social_state +
                              subject +
                              #ehc_inverse_sd_z:social_state +
                              #ehc_inverse_sd_z:subject + 
                              social_state:subject +
                              #ehc_inverse_sd_z:subject:social_state +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)

activity_mod_ESsS <- brm(formula = foraging_initiation ~
                              ehc_inverse_sd_z +
                              social_state +
                              subject +
                              #ehc_inverse_sd_z:social_state +
                              #ehc_inverse_sd_z:subject + 
                              #social_state:subject +
                              #ehc_inverse_sd_z:subject:social_state +
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
                              prior(normal(-2.5, 1.5), class = "Intercept")
                              #prior(exponential(1), class = "sd")
                            ),
                            save_pars = save_pars(all=TRUE),
                            #sample_prior = 'only',
                            family = bernoulli(link = "logit"),
                            control = list(adapt_delta = 0.98)
)



modcompare <- loo_compare(
  loo(activity_mod_ESsSExSsExSSxSExSxSs, k_threshold = 0.7),
  loo(activity_mod_ESsSExSsExSSxS, k_threshold = 0.7),
  loo(activity_mod_ESsSExSsExS, k_threshold = 0.7),
  loo(activity_mod_ESsSExSsSsxS, k_threshold = 0.7),
  loo(activity_mod_ESsSExSSsxS, k_threshold = 0.7),
  loo(activity_mod_ESsSExSs, k_threshold = 0.7),
  loo(activity_mod_ESsSExS, k_threshold = 0.7),
  loo(activity_mod_ESsSSsxS, k_threshold = 0.7),
  loo(activity_mod_ESsS, k_threshold = 0.7))

modcompare <- cbind(modcompare, -2*modcompare[,1], 2*modcompare[,2]) #calc LOOIC difference
colnames(modcompare)[9:10] <- c("looic_diff","se_looic_diff")
print(modcompare[,c("looic","looic_diff","se_looic_diff")], simplify=FALSE, digits=4)
write.csv(as.data.frame(modcompare), "modcompare.csv")
  
summary(activity_mod_ESsSExSsExSSxS)


#POSTERIOR CHECKS
summary(activity_mod_ESsSExSsExSSxS)
pp_check(activity_mod_ESsSExSsExSSxS, type = "bars")
pp_check(activity_mod_ESsSExSsExSSxS, type="dens_overlay",ndraws=1000)

summary(activity_mod)

#plot species and ehc
newd <- data.frame(ehc_inverse_sd_z = rep(seq(min(foraging_initiation_data$ehc_inverse_sd_z), max(foraging_initiation_data$ehc_inverse_sd_z), length.out = 50), 4), subject = factor(rep(c("Surgeonfish","Parrotfish"), each = 100)), social_state = factor(rep(c("Social","Asocial"), each=50, times = 2)))
pmu <- posterior_epred(activity_mod_ESsSExS, newdata = newd, re_formula=NA)

pmu_surgeonfish <- (pmu[, 1:50] + pmu[, 51:100]) / 2
pmu_parrotfish  <- (pmu[, 101:150] + pmu[, 151:200]) / 2

pmu_avg <- cbind(pmu_surgeonfish, pmu_parrotfish)

# Prediction data corresponding to the averaged predictions
pred_newd <- data.frame(
  ehc_inverse_sd_z = rep(
    seq(
      min(foraging_initiation_data$ehc_inverse_sd_z),
      max(foraging_initiation_data$ehc_inverse_sd_z),
      length.out = 50
    ),
    2
  ),
  subject = factor(
    rep(c("Surgeonfish", "Parrotfish"), each = 50)
  )
)

mnmu <- colMeans(pmu_avg)
n <- ncol(pmu_avg)
mean_intervals <- data.frame(mulo95=rep(NA,n), muhi95=rep(NA,n))
for ( i in 1:n ) {
  mean_intervals[i,] <- hpdi(pmu_avg[,i], prob=0.95)
}


preds <- cbind(
  pred_newd,
  mnmu,
  mean_intervals
)

orig_mean <- attr(foraging_initiation_data$ehc_inverse_sd_z, "scaled:center")
orig_sd  <- attr(foraging_initiation_data$ehc_inverse_sd_z, "scaled:scale")

preds$ehc_inverse_sd <- preds$ehc_inverse_sd_z * orig_sd + orig_mean

#the visualisation is where you "undo" the standardisation for the reader by putting the x-axis back in original units
preds |> 
  ggplot() +
  geom_jitter(data = foraging_initiation_data, 
             aes(x = ehc_inverse_sd, y = foraging_initiation, color = subject), 
             alpha = 0.2, size = 1, height = 0.01) +  # fade points back
  geom_ribbon(aes(x = ehc_inverse_sd, ymin = mulo95, ymax = muhi95, fill = subject),  #95% posterior credible intervals (specifically a HPDI), uncertainty around mean 
              alpha = 0.4) +
  geom_line(aes(x = ehc_inverse_sd, y = mnmu, color = subject), 
            linewidth = 1) +
  labs(x = "Experienced Habitat Complexity", 
       y = "Probability initiating a foraging bout") +
  theme_classic(base_size = 14) +
  xlim(0,1)

#plot group size and species
# Prediction data for modifier 2
pred_newd <- data.frame(
  modifier_2_z = rep(
    seq(
      min(foraging_initiation_data$modifier_2_z),
      max(foraging_initiation_data$modifier_2_z),
      length.out = 50
    ),
    2
  ),
  subject = factor(
    rep(c("Surgeonfish", "Parrotfish"), each = 50),
    levels = c("Parrotfish", "Surgeonfish")
  ),
  ehc_inverse_sd_z = 0
)

# Posterior predictions
pmu_avg <- posterior_epred(
  activity_mod,
  newdata = pred_newd,
  re_formula = NA
)

# Posterior mean
mnmu <- colMeans(pmu_avg)

# 95% HPDI
n <- ncol(pmu_avg)

mean_intervals <- data.frame(
  mulo95 = rep(NA, n),
  muhi95 = rep(NA, n)
)

for (i in 1:n) {
  mean_intervals[i, ] <- hpdi(
    pmu_avg[, i],
    prob = 0.95
  )
}

# Combine
preds <- cbind(
  pred_newd,
  mnmu,
  mean_intervals
)

# Undo standardisation
orig_mean <- attr(
  foraging_initiation_data$modifier_2_z,
  "scaled:center"
)

orig_sd <- attr(
  foraging_initiation_data$modifier_2_z,
  "scaled:scale"
)

preds$modifier_2 <- (
  preds$modifier_2_z * orig_sd + orig_mean
)

# Plot
preds |> 
  ggplot() +
  
  geom_jitter(
    data = foraging_initiation_data,
    aes(
      x = modifier_2,
      y = foraging_initiation,
      colour = subject
    ),
    alpha = 0.2,
    size = 1,
    height = 0.01
  ) +
  
  geom_ribbon(
    aes(
      x = modifier_2,
      ymin = mulo95,
      ymax = muhi95,
      fill = subject
    ),
    alpha = 0.4
  ) +
  
  geom_line(
    aes(
      x = modifier_2,
      y = mnmu,
      colour = subject
    ),
    linewidth = 1
  ) +
  
  labs(
    x = "Group size",
    y = "Probability initiating a foraging bout"
  ) +
  
  scale_colour_manual(
    values = c(
      "Parrotfish" = "#E69F00",
      "Surgeonfish" = "#0072B2"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Parrotfish" = "#E69F00",
      "Surgeonfish" = "#0072B2"
    )
  ) +
  
  theme_classic(base_size = 14)


#DISTRIBUTION OF HRI PLOT#DISTRIBUTION OF HRI PLOTsubject
ggplot(
  foraging_initiation_data,
  aes(
    x = factor(foraging_initiation),
    y = ehc_inverse_sd,
    fill = factor(foraging_initiation)
  )
) +
  geom_violin(alpha = 0.3, trim = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.35) +
  scale_fill_manual(
    values = c("0" = "grey60", "1" = "firebrick"),
    guide = "none"
  ) +
  scale_x_discrete(
    labels = c("0" = "No foraging initiation", "1" = "Foraging initiation")
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    colour = "black",
    linewidth = 0.5
  ) +
  labs(
    x = NULL,
    y = "HRI"
  ) +
  theme_classic()

mean(foraging_initiation_data$hri_1)
sd(foraging_initiation_data$hri_1)
# simulate posterior predictive draws
preds <- posterior_predict(activity_mod, ndraws = 1000)  # matrix: draws x observations

sim_resp <- createDHARMa(
  simulatedResponse = t(preds),
  observedResponse  = foraging_initiation_data$foraging_initiation,
  fittedPredictedResponse = colMeans(preds),
  integerResponse = TRUE
)

foraging_initiation_data$resid <- residuals(sim_resp, quantileFunction = qnorm)

acf_by_trial <- foraging_initiation_data %>%
  arrange(ObsID, sec) %>%
  group_by(ObsID) %>%
  group_modify(~{
    a <- acf(.x$resid, plot = FALSE, lag.max = 10, na.action = na.pass)
    tibble(lag = a$lag[,1,1], acf = a$acf[,1,1])
  })

avg_acf <- acf_by_trial %>%
  group_by(lag) %>%
  summarise(mean_acf = mean(acf), se = sd(acf)/sqrt(n())) %>%
  filter(lag > 0)   # lag 0 is always 1, not informative

ggplot(avg_acf, aes(lag, mean_acf)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_acf - 2*se, ymax = mean_acf + 2*se), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Mean residual ACF across trials", x = "Lag (timesteps)")

#max |ACF| = 0.1, no lags exceeded the 95% significance band

#Plotting
newd <- data.frame( hri_z = rep( seq( from = min(foraging_initiation_data$hri_z, na.rm = TRUE), max(foraging_initiation_data$hri_z, na.rm = TRUE), length.out = 50 ), 4 ), social_state = factor( rep(c("Social", "Asocial"), each = 50, times = 2), levels = levels(foraging_initiation_data$social_state) ), subject = factor( rep(c("Parrotfish", "Surgeonfish"), each = 100), levels = levels(foraging_initiation_data$subject) ) )
pmu <- posterior_epred(activity_mod, newdata = newd, re_formula=NA)
ppd <- posterior_predict(activity_mod, newdata = newd, re_formula=NA)
pmu <- (pmu[, 1:50] + pmu[, 51:100]) / 2
mnmu <- colMeans(pmu)
n <- ncol(pmu)
mean_intervals <- data.frame(mulo95=rep(NA,n), muhi95=rep(NA,n))
for ( i in 1:n ) {
  mean_intervals[i,] <- hpdi(pmu[,i], prob=0.95)
}
prediction_intervals <- predictive_interval(ppd, prob=0.95)
prediction_intervals <- data.frame(prediction_intervals)
names(prediction_intervals) <- c("ppdlow95", "ppdhi95")

preds <- cbind(newd, mnmu, mean_intervals, prediction_intervals)

preds |> 
  ggplot() +
  geom_jitter(data = foraging_initiation_data, 
             aes(x = hri_z, y = foraging_initiation), 
             alpha = 0.3, size = 1, height = 0.01, width = 0.01) +  # fade points back
  geom_ribbon(aes(x = hri_z, ymin = mulo95, ymax = muhi95),  #95% posterior credible intervals (specifically a HPDI), uncertainty around mean 
              alpha = 0.4) +
  geom_line(aes(x = hri_z, y = mnmu), 
            linewidth = 1) +
  #geom_line(aes(x = hri, y = ppdlow95, col = social_state), lty = 2) + #prediction intervals, uncertainty aboyt a new observed Y at that X
  #geom_line(aes(x = hri, y = ppdhi95, col = social_state), lty = 2) + #prediction intervals
  labs(x = "HRI", 
       y = "Probability of initiating a foraging bout") +
  theme_bw()


##Model 2 - P(activity state)
stripped_data$activity_state <- as.factor(stripped_data$activity_state)

stripped_data$social_state <- as.factor(stripped_data$social_state)

stripped_data$activity_state_num <- -(as.numeric(stripped_data$activity_state)-2)

stripped_data$ehc_inverse_sd_z <- scale(stripped_data$ehc_inverse_sd)

stripped_data_for_state <- stripped_data |>
  tidyr::drop_na(activity_state_num, ehc_inverse_sd_z, social_state, subject)

activity_mod <- brm(formula = activity_state ~
      ehc_inverse_sd_z +
      social_state +
      #ehc_inverse_sd_z:social_state +
      #ehc_inverse_sd_z:subject +
      #subject:social_state +
      #ar(time = sec, gr = ObsID, p = 1, cov = TRUE) +
        (1|ObsID) +
        (1|site) +
        (1|folder_date),
      data = stripped_data,
      chains = 4,
      cores = 4,
      iter = 4000,
      warmup = 2000,
      prior = c(
        prior(normal(0, 2), class = "b"),
        prior(normal(0, 2), class = "Intercept")
      ),
      save_pars = save_pars(all=TRUE),
      #sample_prior = 'only',
      family = bernoulli(link = "logit"),
      control = list(adapt_delta = 0.98)
)

state_mod <- stan_glmer(activity_state_num ~
             ehc_inverse_sd_z +
             social_state +
             subject +
             ehc_inverse_sd_z:social_state +
             ehc_inverse_sd_z:subject +
             subject:social_state+        
               (1|ObsID) +
               (1|site) +
               (1|folder_date), data = stripped_data, family = binomial(), adapt_delta = 0.95, chains = 4, cores = 4, iter = 3000, seed = 123)
summary(state_mod)

stripped_data %>%
  ggplot(aes(x=ehc_inverse_sd, y = activity_state_num)) +
  geom_jitter(height = 0.01, width = 0.01) +
  geom_smooth()


#PRIOR CHECKS
prior_pred <- posterior_predict(activity_mod)
dim(prior_pred)
prior_rate <- apply(prior_pred, 1, mean)

quantile(
  prior_rate,
  probs = c(.01, .05, .25, .5, .75, .95, .99)
)


#Posterior checks
summary(activity_mod)
plot(activity_mod)
pp_check(state_mod, type="dens_overlay",ndraws=500)
pp_check(state_mod, type = "bars")
pairs(activity_mod)

ci <- posterior_interval(state_mod, prob = 0.95)
overlaps_zero <- ci[,1] < 0 & ci[,2] > 0
data.frame(ci, overlaps_zero)

plot(bout_bites/bout_duration ~ as.factor(subject), data = foraging_bout_stats)



#with interaction
newd <- data.frame(ehc_inverse_sd_z = rep(seq(from = min(stripped_data$ehc_inverse_sd_z, na.rm=TRUE), max(stripped_data$ehc_inverse_sd_z,, na.rm=TRUE), length.out = 50), 2), social_state = factor(rep(c("Asocial", "Social"), each=50)), subject = factor(rep(c("Surgeonfish", "Parrotfish"), each=50)))
pmu <- posterior_epred(activity_mod, newdata = newd, re_formula=NA)
ppd <- posterior_predict(activity_mod, newdata = newd, re_formula=NA)
pmu <- (pmu[, 1:50] + pmu[, 51:100]) / 2
mnmu <- colMeans(pmu)
n <- ncol(pmu)
mean_intervals <- data.frame(mulo95=rep(NA,n), muhi95=rep(NA,n))
for ( i in 1:n ) {
  mean_intervals[i,] <- hpdi(pmu[,i], prob=0.95)
}
prediction_intervals <- predictive_interval(ppd, prob=0.95)
prediction_intervals <- data.frame(prediction_intervals)
names(prediction_intervals) <- c("ppdlow95", "ppdhi95")

preds <- cbind(newd, mnmu, mean_intervals, prediction_intervals)

preds |> 
  ggplot() +
  geom_jitter(data = stripped_data, 
             aes(x = hri_z, y = activity_state_num), 
             alpha = 0.3, size = 1, height = 0.01, width = 0.01) +  # fade points back
  geom_ribbon(aes(x = hri_z, ymin = mulo95, ymax = muhi95),  #95% posterior credible intervals (specifically a HPDI), uncertainty around mean 
              alpha = 0.4) +
  geom_line(aes(x = hri_z, y = mnmu), 
            linewidth = 1) +
  #geom_line(aes(x = hri, y = ppdlow95, col = social_state), lty = 2) + #prediction intervals, uncertainty aboyt a new observed Y at that X
  #geom_line(aes(x = hri, y = ppdhi95, col = social_state), lty = 2) + #prediction intervals
  labs(x = "HRI", 
       y = "Probability of being in foraging state") +
  theme_bw()


#withotu interaction
newd <- data.frame(hri_z = rep(seq(from = min(stripped_data$hri_z, na.rm=TRUE), max(stripped_data$hri_z, na.rm=TRUE), length.out = 100)), sec   = 1, ObsID = paste0("new_obs_", 1:100) )
pmu <- posterior_epred(activity_mod, newdata = newd, re_formula=NA, allow_new_levels = TRUE)
mnmu <- colMeans(pmu)
n <- ncol(pmu)
mean_intervals <- data.frame(mulo95=rep(NA,n), muhi95=rep(NA,n))
for ( i in 1:n ) {
  mean_intervals[i,] <- hpdi(pmu[,i], prob=0.95)
}


preds <- cbind(newd, mnmu, mean_intervals)

preds |> 
  ggplot() +
  geom_point(data = stripped_data, 
             aes(x = hri_z, y = activity_state_num), 
             alpha = 0.3, size = 1) +  # fade points back
  geom_ribbon(aes(x = hri_z, ymin = mulo95, ymax = muhi95),  #95% posterior credible intervals (specifically a HPDI), uncertainty around mean 
              alpha = 0.4) +
  geom_line(aes(x = hri_z, y = mnmu), 
            linewidth = 1) +
  #geom_line(aes(x = hri, y = ppdlow95, col = social_state), lty = 2) + #prediction intervals, uncertainty aboyt a new observed Y at that X
  #geom_line(aes(x = hri, y = ppdhi95, col = social_state), lty = 2) + #prediction intervals
  labs(x = "HRI", 
       y = "Activity state") +
  theme_bw()

stripped_data$Bites_sc <- scale(stripped_data$Bites)

##Model 3 - Bite count
bite_mod <- brm(formula = Bites ~
                      hri_z +
                      social_state +
                      hri_z:social_state +
                      subject + 
                      #ar(time = sec, gr = ObsID, p = 1, cov = TRUE) +
                      (1|ObsID) +
                      (1|site) +
                      (1|folder_date),
                    data = stripped_data,
                    chains = 4,
                    cores = 4,
                    iter = 4000,
                    warmup = 2000,
                    prior = c(
                      prior(normal(0, 2), class = "b"),
                      prior(normal(0, 2), class = "Intercept")
                    ),
                    save_pars = save_pars(all=TRUE),
                    sample_prior = 'only',
                    family = poisson(),
                    control = list(adapt_delta = 0.98)
)

#is var > mean, then overdispersed and poisson violated, use nb
mean(stripped_data$Bites)
var(stripped_data$Bites)
summary(bite_mod)
plot(bite_mod)
pp_check(bite_mod, type="dens_overlay",ndraws=500)
pp_check(bite_mod, type = "bars")

#plotting bites
newd <- data.frame(hri_z = rep(seq(from = min(stripped_data$hri_z, na.rm=TRUE), max(stripped_data$hri_z, na.rm=TRUE), length.out = 50), 4), 
                   social_state = factor(rep(c("Asocial", "Social"), each=100)), 
                   subject = factor(rep(c("Surgeonfish", "Parrotfish"), each=50)))
pmu <- posterior_epred(bites_mod, newdata = newd, re_formula=NA)
ppd <- posterior_predict(bites_mod, newdata = newd, re_formula=NA)

# average over subject 
pmu_asocial <- (pmu[, 1:50] + pmu[, 51:100]) / 2 
pmu_social <- (pmu[, 101:150] + pmu[, 151:200]) / 2 
mnmu <- c(colMeans(pmu_asocial), colMeans(pmu_social)) 
n <- ncol(pmu) 
mean_intervals <- data.frame(mulo95=rep(NA,n), muhi95=rep(NA,n)) 

for (i in 1:n) { 
  mean_intervals[i,] <- hpdi(pmu[,i], prob=0.95) 
  } 
n <- ncol(pmu_social) 
mean_intervals_social <- data.frame(mulo95=rep(NA,n), muhi95=rep(NA,n)) 
for (i in 1:n) { 
  mean_intervals_social[i,] <- hpdi(pmu_social[,i], prob=0.95) 
} 

prediction_intervals <- predictive_interval(ppd, prob=0.95) 
prediction_intervals <- data.frame(prediction_intervals) names(prediction_intervals) <- c("ppdlow95", "ppdhi95")

pmu <- (pmu[, 1:50] + pmu[, 51:100]) / 2
mnmu <- colMeans(pmu)
n <- ncol(pmu)
mean_intervals <- data.frame(mulo95=rep(NA,n), muhi95=rep(NA,n))
for ( i in 1:n ) {
  mean_intervals[i,] <- hpdi(pmu[,i], prob=0.95)
}
prediction_intervals <- predictive_interval(ppd, prob=0.95)
prediction_intervals <- data.frame(prediction_intervals)
names(prediction_intervals) <- c("ppdlow95", "ppdhi95")

preds <- cbind(newd, mnmu, mean_intervals, prediction_intervals)

preds |> 
  ggplot() +
  geom_jitter(data = stripped_data, 
              aes(x = hri_z, y = activity_state_num), 
              alpha = 0.3, size = 1, height = 0.01, width = 0.01) +  # fade points back
  geom_ribbon(aes(x = hri_z, ymin = mulo95, ymax = muhi95),  #95% posterior credible intervals (specifically a HPDI), uncertainty around mean 
              alpha = 0.4) +
  geom_line(aes(x = hri_z, y = mnmu), 
            linewidth = 1) +
  #geom_line(aes(x = hri, y = ppdlow95, col = social_state), lty = 2) + #prediction intervals, uncertainty aboyt a new observed Y at that X
  #geom_line(aes(x = hri, y = ppdhi95, col = social_state), lty = 2) + #prediction intervals
  labs(x = "HRI", 
       y = "Probability of being in foraging state") +
  theme_bw()


#Model 4 - Bite rate
#checking if nbinom is justified
mean(bout_stats$bout_bites)
var(bout_stats$bout_bites)

bout_stats <- bout_stats |>
  tidyr::drop_na(bout_bites, bout_duration, bout_ehc_inverse_sd_z, social_state, bout_type)

foraging_bout_stats <- bout_stats %>%
  filter(bout_type == "Foraging")

bite_rate_mod <- brm(formula = bout_bites ~
                       bout_ehc_inverse_sd_z +
                      social_state +
                      subject +
                      bout_ehc_inverse_sd_z:social_state +
                      bout_ehc_inverse_sd_z:subject + 
                      social_state:subject +
                      offset(log(bout_duration)) +
                      #ar(time = sec, gr = ObsID, p = 1, cov = TRUE) +
                      (1|Trial_ID) +
                      (1|site) +
                      (1|folder_date),
                    data = foraging_bout_stats,
                    chains = 4,
                    cores = 4,
                    iter = 4000,
                    warmup = 2000,
                    prior = c(
                      prior(normal(0, 2), class = "b"),
                      prior(normal(0, 2), class = "Intercept")
                    ),
                    save_pars = save_pars(all=TRUE),
                    sample_prior = 'only',
                    family = negbinomial(),
                    control = list(adapt_delta = 0.98)
)

bout_bite_model <- stan_glmer.nb(bout_bites ~
             bout_ehc_inverse_sd_z +
             social_state +
             subject +
             bout_ehc_inverse_sd_z:social_state +
             bout_ehc_inverse_sd_z:subject +
             social_state:subject +
               social_state:subject:bout_ehc_inverse_sd_z +  
             offset(log(bout_duration)) +
             #ar(time = sec, gr = ObsID, p = 1, cov = TRUE) +
             (1|Trial_ID) +
             (1|site) +
             (1|folder_date),
           data = foraging_bout_stats, adapt_delta = 0.99)

print(summary(bout_bite_model)[,c("mean","sd","n_eff","Rhat")], digits=3)

launch_shinystan(bout_bite_model)
plot(bout_bite_model, "hist", bins=50)
pp_check(bout_bite_model, plotfun = "dens_overlay")

pp_check(bout_bite_model, type="dens_overlay",ndraws=500)
posterior_interval(bout_bite_model, prob = 0.95)

ci <- posterior_interval(bout_bite_model, prob = 0.95)
overlaps_zero <- ci[,1] < 0 & ci[,2] > 0
data.frame(ci, overlaps_zero)

plot(bout_bites/bout_duration ~ as.factor(subject), data = foraging_bout_stats)


#Distribution
#Heights
hri_data_f_height_f <- full_data %>%
  select(ObsID, hit_rate, miss_rate, frame.x, subject, starts_with("height")) %>%
  mutate(hit_rate = as.numeric(hit_rate),
          across(height_bin_0_0_0_2m:height_bin_3_6_3_8m, ~ as.numeric(as.character(.x))),
         sum_bins = rowSums(across(height_bin_0_0_0_2m:height_bin_3_6_3_8m), na.rm=TRUE),
         across(height_bin_0_0_0_2m:height_bin_3_6_3_8m, ~.x * (hit_rate/100) * 50000, .names = "ct_{.col}"))


hri_data_f_height_f_long <- hri_data_f_height_f %>%
  select(!c(height_bin_0_0_0_2m:height_bin_3_6_3_8m, sum_bins)) %>%
  pivot_longer(cols = c(ct_height_bin_0_0_0_2m:ct_height_bin_3_6_3_8m),
               names_to = "Bins", values_to = "Counts", names_prefix = "ct_distance_bin_") %>%
  mutate(Bins_num = as.numeric(as.factor(Bins))) %>%
  mutate(Counts = round(Counts, 0),
         Proportions = Counts / ((hit_rate/100)*50000)) %>%
  mutate(Bins_center = case_when(
             Bins_num == 1 ~ 0.1,
             Bins_num == 2 ~ 0.3,
             Bins_num == 3 ~ 0.5,
             Bins_num == 4 ~ 0.7,
             Bins_num == 5 ~ 0.9,
             Bins_num == 6 ~ 1.1,
             Bins_num == 7 ~ 1.3,
             Bins_num == 8 ~ 1.4,
             Bins_num == 9 ~ 1.5,
             Bins_num == 10 ~ 1.7,
             Bins_num == 11 ~ 1.9,
             Bins_num == 12 ~ 2.1,
             Bins_num == 13 ~ 2.3,
             Bins_num == 14 ~ 2.5,
             Bins_num == 15 ~ 2.7,
             Bins_num == 16 ~ 2.9,
             Bins_num == 17 ~ 3.1,
             Bins_num == 18 ~ 3.3,
             Bins_num == 19 ~ 3.5,
           )) %>%
  mutate(Bins_num = as.numeric(Bins_num))

#Distances
hri_data_f_distance_f <- full_data %>%
  select(ObsID, hit_rate, miss_rate, frame.x, subject, starts_with("distance")) %>%
  mutate(hit_rate = as.numeric(hit_rate),
         across(distance_bin_0_0_0_5m:distance_bin_1_5_2_0m, ~ as.numeric(as.character(.x))),
         sum_bins = rowSums(across(distance_bin_0_0_0_5m:distance_bin_1_5_2_0m), na.rm=TRUE),
         across(distance_bin_0_0_0_5m:distance_bin_1_5_2_0m, ~.x * (hit_rate/100) * 50000, .names = "ct_{.col}"))


hri_data_f_distance_f_long <- hri_data_f_distance_f %>%
  select(!c(distance_bin_0_0_0_5m:distance_bin_1_5_2_0m, sum_bins)) %>%
  pivot_longer(cols = c(ct_distance_bin_0_0_0_5m:ct_distance_bin_1_5_2_0m),
               names_to = "Bins", values_to = "Counts", names_prefix = "ct_distance_bin_") %>%
  mutate(Bins_num = as.numeric(as.factor(Bins))) %>%
  mutate(Counts = round(Counts, 0),
         Proportions = Counts / ((hit_rate/100)*50000)) %>%
  mutate(Bins_center = case_when(
    Bins_num == 1 ~ 0.25,
    Bins_num == 2 ~ 0.75,
    Bins_num == 3 ~ 1.25,
    Bins_num == 4 ~ 1.75
  )) %>%
  mutate(Bins_num = as.numeric(Bins_num))


#Hit rate
per_frame <- hri_data_f_distance_f_long %>%
  group_by(ObsID, frame.x) %>%
  reframe(hit_rate = first(hit_rate), subject = first(subject)) %>%
  mutate(frame.x = as.numeric(frame.x), ObsID = as.factor(ObsID))

#hit rate densities per trial
dist.hits <- ggplot(data = per_frame) +
  geom_density(aes(x = hit_rate, color = as.factor(subject), group = ObsID), fill = NA) +
  labs(x = "Hit Rate", y = "Density") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") 
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/hit_rate.png", plot = dist.hits, width = 10, height = 5)

#hit rate across the trial
ggplot(per_frame, aes(x = frame.x, y = hit_rate, color = as.factor(ObsID))) +
  geom_line(aes(group = ObsID, alpha = 0.1)) +
  labs(x = "Frame", y = "Hit Rate") +
  theme_classic(base_size = 14)

dist.bins <- ggplot(data = hri_data_f_distance_f_long) +
  geom_boxplot(aes(x = Bins, y = Counts), fill = NA) +
  geom_point(aes(x = Bins, y = Counts, color = as.factor(ObsID))) +
  theme_classic()
dist.bins

#Hit heights
height.dens <- ggplot(data = hri_data_f_height_f_long %>% filter(!is.na(Counts))) +
  geom_density(aes(x = Bins_center, weight = Counts, color = as.factor(subject), group = ObsID), fill = NA) +
  # scale_x_continuous(breaks = seq(0.25, 2, 0.5)) +
  scale_y_continuous(breaks = seq(0, 6, 0.1)) +
  labs(x = "Hit heights (m)", y = "Density") +
  xlim(c(0,2)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") 
height.dens
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/height.dens.png", plot = height.dens, width = 10, height = 5)


#Hit distances
dist.dens <- ggplot(data = hri_data_f_distance_f_long %>% filter(!is.na(Counts))) +
  geom_density(aes(x = Bins_center, weight = Counts, color = as.factor(subject), group = ObsID), fill = NA, bw = 0.3) +
  # scale_x_continuous(breaks = seq(0.25, 2, 0.5)) +
  scale_y_continuous(breaks = seq(0, 6, 0.1)) +
  labs(x = "Hit distances (m)", y = "Density") +
  #xlim(c(0.25,2)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") 
dist.dens
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/dist.dens.png", plot = dist.dens, width = 10, height = 5)

#EHC
ehc <- full_data %>%
  ggplot(aes(x = ehc, group = as.factor(ObsID), color = as.factor(subject))) + 
  geom_density() +
  #xlim(0,5) +
  labs(x = "Experienced Habitat Complexity", y = "Density") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") 

ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/ehc.png", plot = ehc, width = 10, height = 5)

  
ehc_inverse_sd <- full_data %>%
  ggplot(aes(x = ehc_inverse_sd, group = as.factor(ObsID), color = as.factor(subject))) + 
  geom_density() +
  #xlim(0,5) +
  labs(x = "Experienced Habitat Complexity", y = "Density") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") 
ggsave("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/paper/figures/ehc_inverse_sd.png", plot = ehc_inverse_sd, width = 10, height = 5)


