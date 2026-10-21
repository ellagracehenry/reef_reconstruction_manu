levels(full_data_strip$subject)
library(posterior)
library(bayestestR)

mod <- activity_models$lESSs_lExS_lExSs_SxSs_lExSxSs   # replace with your model object
d   <- as_draws_df(mod)

main <- d$`b_log_ehc_z:social_stateSocial`
int3 <- d$`b_log_ehc_z:subjectSurgeonfish:social_stateSocial`

# If levels() shows Surgeonfish is the SECOND level, swap the signs:
# sf_social_minus_solitary    <- main + s * int3
# other_social_minus_solitary <- main - s * int3

contrasts <- data.frame(
  refernece_species = main,
  surgeonfish = main + int3,
  species_differences = int3
)

describe_posterior(
  contrasts,
  ci = 0.95,
  test = c("pd", "rope"),
  rope_range = c(-0.1, 0.1)
)

d <- as_draws_df(activity_models$lESSs_lExS_lExSs_SxSs_lExSxSs)

sf_social_minus_solitary <- d$`b_log_ehc_z:social_stateSocial` +
  d$`b_log_ehc_z:subjectSurgeonfish:social_stateSocial`

plot(sf_social_minus_solitary)
