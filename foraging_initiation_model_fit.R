library(brms)
library(loo)
library(dplyr)
library(purrr)
library(stringr)


# ============================================================
# VARIABLES
# ============================================================

vars <- c(
  E = "ehc_z",
  G = "modifier_2_z",
  M = "modifier_1",
  S = "subject"
)


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
        "ehc_z",
        "E"
      )
      
      x <- str_replace_all(
        x,
        "group_size_binned",
        "G"
      )
      
      x <- str_replace_all(
        x,
        "modifier_1",
        "M"
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
      "foraging_initiation ~",
      paste(
        terms,
        collapse = " + "
      ),
      "+ (1|ObsID) + (1|site) + (1|folder_date)"
    )
  )
  
  
  brm(
    formula = formula,
    
    data = foraging_initiation_data,
    
    chains = 4,
    cores = 4,
    
    iter = 4000,
    warmup = 2000,
    
    prior = c(
      
      prior(
        normal(0, 1),
        class = "b"
      ),
      
      prior(
        normal(-2.5, 1.5),
        class = "Intercept"
      )
      
      # prior(
      #   exponential(1),
      #   class = "sd"
      # )
      
    ),
    
    save_pars = save_pars(
      all = TRUE
    ),
    
    family = bernoulli(
      link = "logit"
    ),
    
    control = list(
      adapt_delta = 0.98
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
  
}


# ============================================================
# SAVE ALL MODELS
# ============================================================

save(
  activity_models,
  model_terms,
  model_names,
  additive_models,
  interaction_models,
  file = "activity_models_all.RData"
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
  activity_models$E_add
)
