
library(brms)
library(loo)
library(dplyr)
library(purrr)
library(stringr)

#sort factors
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

full_data_strip <- full_data |>
  tidyr::drop_na(ehc_z, subject, modifier_1, modifier_2_z, bites_binary) 

plot(bites_binary ~ log_ehc_z_squared, full_data_strip )

class(full_data_strip$Trial_ID)
hist(full_data$log_ehc_z)
#Test
activity_mod_ehc_filtered <- brm(
  formula = bites_binary ~
    ehc_z*subject*social_state +
   (1 | Trial_ID) +
    (1 | site) +
    (1 | folder_date),
  
  data = full_data_strip,
  
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
  
  #sample_prior = "only",
  
  family = bernoulli(link = "logit"),
  
  control = list(
    adapt_delta = 0.99
  )
)


plot(activity_mod)
pp_check(activity_mod_log_ehc_filtered, ndraws=100)
pairs(activity_mod)
summary(activity_mod_log_ehc_filtered)
prior_summary(activity_mod1)
# ============================================================
# VARIABLES
# ============================================================

vars <- c(
  lE = "log_ehc_z",
 # lE2 = "log_ehc_z_squared",
  S = "subject",
  Ss = "social_state"
)

class(full_data_strip$log_ehc_z)
hist(full_data_strip$log_ehc_z)
hist(full_data_strip$log_ehc_z_squared)
class(full_data_strip$social_state)
# ============================================================
# ALL FACTOR SUBSETS
# ============================================================

# Generate every non-empty subset of E, G, M, S
factor_sets <- unlist(
  lapply(
    1:length(vars),
    function(k) {
      combn(
        names(vars),
        k,
        simplify = FALSE
      )
    }
  ),
  recursive = FALSE
)


# ============================================================
# FUNCTION: CONVERT FACTOR NAMES TO VARIABLES
# ============================================================

factor_vars <- function(factors) {
  unname(vars[factors])
}


# ============================================================
# FUNCTION: CREATE MAIN-EFFECTS-ONLY MODEL
# ============================================================

make_additive_terms <- function(factors) {
  
  factor_vars(factors)
  
}


# ============================================================
# FUNCTION: CREATE ALL POSSIBLE INTERACTIONS
# ============================================================

make_all_interactions <- function(factors) {
  
  if (length(factors) < 2) {
    return(character(0))
  }
  
  interactions <- character(0)
  
  for (k in 2:length(factors)) {
    
    interaction_sets <- combn(
      factors,
      k,
      simplify = FALSE
    )
    
    interaction_terms <- vapply(
      interaction_sets,
      function(x) {
        
        paste(
          factor_vars(x),
          collapse = ":"
        )
        
      },
      character(1)
    )
    
    interactions <- c(
      interactions,
      interaction_terms
    )
  }
  
  interactions
  
}


# ============================================================
# FUNCTION: CHECK HIERARCHY
# ============================================================

is_hierarchical <- function(selected_interactions) {
  
  # No interactions = valid model
  if (length(selected_interactions) == 0) {
    return(TRUE)
  }
  
  
  for (term in selected_interactions) {
    
    variables <- strsplit(
      term,
      ":",
      fixed = TRUE
    )[[1]]
    
    k <- length(variables)
    
    
    # Pairwise interactions are always hierarchical
    if (k < 3) {
      next
    }
    
    
    # For higher-order interactions,
    # every lower-order interaction must exist
    for (j in 2:(k - 1)) {
      
      lower_combinations <- combn(
        variables,
        j,
        simplify = FALSE
      )
      
      
      for (lower in lower_combinations) {
        
        lower_term <- paste(
          lower,
          collapse = ":"
        )
        
        
        if (!(lower_term %in% selected_interactions)) {
          return(FALSE)
        }
        
      }
    }
  }
  
  TRUE
  
}


# ============================================================
# FUNCTION: CREATE MODEL NAME
# ============================================================

make_model_name <- function(factors, interactions) {
  
  # ----------------------------------------------------------
  # Factor abbreviation
  # ----------------------------------------------------------
  
  factor_name <- paste(
    factors,
    collapse = ""
  )
  
  
  # ----------------------------------------------------------
  # Additive model
  # ----------------------------------------------------------
  
  if (length(interactions) == 0) {
    
    return(
      paste0(
        factor_name,
        "_add"
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Convert interactions to readable names
  # ----------------------------------------------------------
  
  readable_interactions <- vapply(
    interactions,
    function(x) {
      
      x <- str_replace_all(
        x,
        "log_ehc_z",
        "lE"
      )
      
     # x <- str_replace_all(
    #    x,
    #    "log_ehc_z_squared",
    #    "lE2"
    #  )
      
      x <- str_replace_all(
        x,
        "social_state",
        "Ss"
      )
      
      x <- str_replace_all(
        x,
        "subject",
        "S"
      )
      
      x <- str_replace_all(
        x,
        ":",
        "x"
      )
      
      x
      
    },
    character(1)
  )
  
  
  paste0(
    factor_name,
    "_",
    paste(
      readable_interactions,
      collapse = "_"
    )
  )
  
}


# ============================================================
# CREATE ALL ADDITIVE MODELS
# ============================================================

additive_models <- list()


for (factors in factor_sets) {
  
  terms <- make_additive_terms(factors)
  
  name <- make_model_name(
    factors,
    character(0)
  )
  
  additive_models[[name]] <- terms
  
}


# ============================================================
# CREATE ALL HIERARCHICAL INTERACTION MODELS
# ============================================================

interaction_models <- list()


for (factors in factor_sets) {
  
  # ----------------------------------------------------------
  # Get every possible interaction for this factor set
  # ----------------------------------------------------------
  
  possible_interactions <- make_all_interactions(
    factors
  )
  
  
  # No interactions possible for 1-factor models
  if (length(possible_interactions) == 0) {
    next
  }
  
  
  n_interactions <- length(
    possible_interactions
  )
  
  
  # ----------------------------------------------------------
  # Generate every subset of possible interactions
  # ----------------------------------------------------------
  
  for (i in 1:n_interactions) {
    
    interaction_subsets <- combn(
      possible_interactions,
      i,
      simplify = FALSE
    )
    
    
    for (selected_interactions in interaction_subsets) {
      
      # ------------------------------------------------------
      # Keep only hierarchical models
      # ------------------------------------------------------
      
      if (!is_hierarchical(
        selected_interactions
      )) {
        next
      }
      
      
      # ------------------------------------------------------
      # Always include all main effects
      # ------------------------------------------------------
      
      terms <- c(
        factor_vars(factors),
        selected_interactions
      )
      
      
      # ------------------------------------------------------
      # Create model name
      # ------------------------------------------------------
      
      name <- make_model_name(
        factors,
        selected_interactions
      )
      
      
      interaction_models[[name]] <- terms
      
    }
  }
}


# ============================================================
# COMBINE ADDITIVE + INTERACTION MODELS
# ============================================================

model_terms <- c(
  additive_models,
  interaction_models
)

model_terms <- model_terms[
  vapply(
    model_terms,
    function(x) {
      !("log_ehc_z_squared" %in% x) || ("log_ehc_z" %in% x)
    },
    logical(1)
  )
]

model_terms <- model_terms[
  vapply(
    model_terms,
    function(x) {
      
      # squared term requires linear term
      squared_hierarchy <- 
        !("log_ehc_z_squared" %in% x) ||
        ("log_ehc_z" %in% x)
      
      # don't allow log_ehc_z × log_ehc_z_squared
      no_cubic <- !any(
        grepl(
          "log_ehc_z:log_ehc_z_squared|log_ehc_z_squared:log_ehc_z",
          x
        )
      )
      
      squared_hierarchy && no_cubic
    },
    logical(1)
  )
]
# ============================================================
# REMOVE ANY DUPLICATE MODELS
# ============================================================

model_signature <- vapply(
  model_terms,
  function(x) {
    
    paste(
      sort(x),
      collapse = " + "
    )
    
  },
  character(1)
)


keep <- !duplicated(
  model_signature
)


model_terms <- model_terms[keep]


# ============================================================
# MODEL NAMES
# ============================================================

model_names <- names(
  model_terms
)


# ============================================================
# CHECK NUMBER OF MODELS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Number of additive models:",
  length(additive_models),
  "\n"
)

cat(
  "Number of interaction models:",
  length(interaction_models),
  "\n"
)

cat(
  "TOTAL UNIQUE MODELS:",
  length(model_terms),
  "\n"
)

cat(
  "============================================\n"
)


# ============================================================
# DISPLAY ALL MODELS
# ============================================================

for (i in seq_along(model_terms)) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    i,
    "of",
    length(model_terms),
    "\n"
  )
  
  cat(
    "Model:",
    model_names[i],
    "\n"
  )
  
  cat(
    "Terms:",
    paste(
      model_terms[[i]],
      collapse = " + "
    ),
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
}


# ============================================================
# BRMS MODEL FUNCTION
# ============================================================

fit_activity_model <- function(terms) {
  
  formula <- as.formula(
    paste(
      "bites_binary ~",
      paste(
        terms,
        collapse = " + "
      ),
      "+ (1|Trial_ID) + (1|site) + (1|folder_date)"
    )
  )
  
  
  brm(
    formula = formula,
    
    data = full_data_strip,
    
    chains = 4,
    cores = 4,
    
    iter = 2000,
    warmup =1000,
    
    prior = c(
      
      prior(
        normal(0, 1),
        class = "b"
      ),
      
      prior(
        normal(-2.5, 1.5),
        class = "Intercept"
      ),
      
       prior(
         exponential(1),
         class = "sd"
       )
      
    ),
    
    save_pars = save_pars(
      all = TRUE
    ),
    
    family = bernoulli(
      link = "logit"
    ),
    
    control = list(
      adapt_delta = 0.995
    )
  )
}


# ============================================================
# FIT ALL MODELS
# ============================================================

activity_models <- vector(
  "list",
  length(model_terms)
)

names(activity_models) <- model_names


for (i in seq_along(model_terms)) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Fitting model",
    i,
    "of",
    length(model_terms),
    "\n"
  )
  
  cat(
    "Model:",
    model_names[i],
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  activity_models[[i]] <- fit_activity_model(
    model_terms[[i]]
  )
  
  save(
    activity_models,
    model_terms,
    model_names,
    additive_models,
    interaction_models,
    file = "bite_models_log_ehc_filtered.RData"
  )
  
}


# ============================================================
# SAVE ALL MODELS
# ============================================================

save(
  bite_models_ehc,
  model_terms,
  model_names,
  additive_models,
  interaction_models,
  file = "bite_models_log_ehc_filtered.RData"
)


# ============================================================
# CALCULATE LOO
# ============================================================

loo_mods <- lapply(
  activity_models,
  function(mod) {
    
    loo(
      mod,
      k_threshold = 0.7
    )
    
  }
)


# ============================================================
# COMPARE MODELS
# ============================================================

modcompare <- loo_compare(
  loo_mods
)


# ============================================================
# PRINT MODEL COMPARISON
# ============================================================

print(
  modcompare
)


# ============================================================
# EXAMPLE: INSPECT A MODEL
# ============================================================

summary(
  activity_models$lESSs_lExS_lExSs_SxSs_lExSxSs 
)

EGMS_sd_zxmodifier_2_z_sd_zxM_sd_zxS_modifier_2_zxS_sd_zxmodifier_2_zxS <- activity_models$EGMS_sd_zxmodifier_2_z_sd_zxM_sd_zxS_modifier_2_zxS_sd_zxmodifier_2_zxS

load("/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/BORIS/frame_idx/activity_models.RData")
activity_models_sd_z <- activity_models

#Figure
# Three way Prediction data corresponding to the averaged predictions
pred_newd <- expand.grid(
  log_ehc_z = seq(
      min(full_data_strip$log_ehc_z),
      max(full_data_strip$log_ehc_z),
      length.out = 50
    ),
  subject = factor(
    c("Surgeonfish", "Parrotfish"), 
    levels = levels(full_data_strip$subject)),
  social_state = factor(
    c("Asocial", "Social"), 
    levels = levels(full_data_strip$social_state))
)

pmu <- posterior_epred(activity_models$lESSs_lExS_lExSs_SxSs, newdata = pred_newd, re_formula=NA)

preds <- cbind(
  pred_newd,
  mnmu = apply(pmu, 2, mean),
  t(apply(pmu, 2, hpdi, prob = 0.95))
)

names(preds)[names(preds) == "1"] <- "mulo95"
names(preds)[names(preds) == "2"] <- "muhi95"

orig_mean <- attr(full_data_strip$ehc_z, "scaled:center")
orig_sd  <- attr(full_data_strip$ehc_z, "scaled:scale")

preds$log_ehc <- preds$log_ehc_z * orig_sd + orig_mean

#the visualisation is where you "undo" the standardisation for the reader by putting the x-axis back in original units
preds |> 
  ggplot() +
  geom_jitter(data = full_data_strip, 
              aes(x = log_ehc, y = bites_binary, color = social_state), 
              alpha = 0.1, size = 1, height = 0.01) +  # fade points back
  geom_ribbon(data = preds, aes(x = log_ehc, ymin = mulo95, ymax = muhi95, fill = social_state),  #95% posterior credible intervals (specifically a HPDI), uncertainty around mean 
              alpha = 0.2) +
  geom_line(data = preds, aes(x = log_ehc, y = mnmu, color = social_state), 
            linewidth = 1) +
  facet_wrap(~subject) +
  labs(x = "Experienced Habitat Complexity", 
       y = "Probability of biting") +
  theme_classic(base_size = 24) #+
  #coord_cartesian(xlim = c(-5,5))

#trying log scale predictor
full_data_strip$log_ehc_z <- scale(log(full_data_strip$ehc))
full_data_strip$log_ehc <- log(full_data_strip$ehc)
hist(full_data_strip$log_ehc_z)

full_data_strip$subject <- factor(
  full_data_strip$subject,
  levels = c("Surgeonfish", "Parrotfish"),
  ordered = FALSE   # makes sure you get treatment coding, not .L
)

levels(full_data_strip$subject)

activity_mod_logehc_relvel_nonf <- brm(
  formula = bites_binary ~
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
    adapt_delta = 0.995
  )
)

save(activity_mod_logehc_relvel_nonf, file = "bites_log_ehc_relevel_nonf.RData")
save(activity_mod_logehc_relvel, file = "bites_log_ehc_relevel_filt.RData")

mod_sf <- update(
  activity_models$lESSs_lExS_lExSs_SxSs_lExSxSs,
  newdata = full_data_strip,
  recompile = FALSE,
  iter = 4000, warmup = 2000, chains = 4, cores = 4
)


summary(activity_mod_logehc_relvel)
summary(activity_models$lESSs_lExS_lExSs_SxSs_lExSxSs)
summary(activity_models$ESSs_ExS_ExSs_SxSs_ExSxSs)
loo(activity_mod_ehc_filtered, activity_mod_log_ehc_filtered)

pp_check(activity_mod_logehc)
pp_check(activity_models$ESSs_ExS_ExSs_SxSs_ExSxSs)

max(full_data_strip$ehc[full_data_strip$social_state == "Social" & full_data_strip$bites_binary == 1])
max(full_data_strip$ehc[full_data_strip$social_state == "Asocial"& full_data_strip$bites_binary == 1])
min(full_data_strip$ehc[full_data_strip$social_state == "Social"& full_data_strip$bites_binary == 1])
min(full_data_strip$ehc[full_data_strip$social_state == "Asocial"& full_data_strip$bites_binary == 1])

#Two way figure
newd <- data.frame(log_ehc_z = rep(seq(from = min(full_data_strip$log_ehc_z, na.rm=TRUE), max(full_data_strip$log_ehc_z, na.rm=TRUE), length.out = 50), 4), 
                   social_state = factor(rep(c("Asocial", "Social"), each=100)), 
                   subject = factor(rep(c("Surgeonfish", "Parrotfish"), each=50)))
pmu <- posterior_epred(activity_models$lESSs_lExS_lExSs_SxSs, newdata = newd, re_formula=NA)
ppd <- posterior_predict(activity_models$lESSs_lExS_lExSs_SxSs, newdata = newd, re_formula=NA)

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
  geom_jitter(data = full_data_strip, 
              aes(x = log_ehc_z, y = bites_binary), 
              alpha = 0.3, size = 1, height = 0.01, width = 0.01) +  # fade points back
  geom_ribbon(aes(x = log_ehc_z, ymin = mulo95, ymax = muhi95),  #95% posterior credible intervals (specifically a HPDI), uncertainty around mean 
              alpha = 0.4) +
  geom_line(aes(x = log_ehc_z, y = mnmu), 
            linewidth = 1) +
  #geom_line(aes(x = hri, y = ppdlow95, col = social_state), lty = 2) + #prediction intervals, uncertainty aboyt a new observed Y at that X
  #geom_line(aes(x = hri, y = ppdhi95, col = social_state), lty = 2) + #prediction intervals
  labs(x = "HRI", 
       y = "Probability of being in foraging state") +
  facet_wrap(~subject)+
  theme_bw()
