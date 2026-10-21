#posterior draws
library(posterior)

draws <- as_draws_df(
  seafloor_bites_model_conspecifics)
draws$ehc_ref_type0 <-
  draws$b_log_ehc_z

draws$ehc_ref_type1 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:social_type1`

draws$ehc_ref_type2 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:social_type2`

draws$ehc_surgeon_type0 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:subjectSurgeonfish`

draws$ehc_surgeon_type1 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:subjectSurgeonfish` +
  draws$`b_log_ehc_z:social_type1` +
  draws$`b_log_ehc_z:subjectSurgeonfish:social_type1`

draws$ehc_surgeon_type2 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:subjectSurgeonfish` +
  draws$`b_log_ehc_z:social_type2` +
  draws$`b_log_ehc_z:subjectSurgeonfish:social_type2`
summarize_posterior <- function(x) {
  c(
    Estimate = median(x),
    Lower_95 = quantile(x, 0.025),
    Upper_95 = quantile(x, 0.975),
    P_positive = mean(x > 0)
  )
}

rbind(
  ref_type0      = summarize_posterior(draws$ehc_ref_type0),
  ref_type1      = summarize_posterior(draws$ehc_ref_type1),
  ref_type2      = summarize_posterior(draws$ehc_ref_type2),
  surgeon_type0  = summarize_posterior(draws$ehc_surgeon_type0),
  surgeon_type1  = summarize_posterior(draws$ehc_surgeon_type1),
  surgeon_type2  = summarize_posterior(draws$ehc_surgeon_type2)
)

rbind(
  ref_type0      = summarize_posterior(exp(draws$ehc_ref_type0)),
  ref_type1      = summarize_posterior(exp(draws$ehc_ref_type1)),
  ref_type2      = summarize_posterior(exp(draws$ehc_ref_type2)),
  surgeon_type0  = summarize_posterior(exp(draws$ehc_surgeon_type0)),
  surgeon_type1  = summarize_posterior(exp(draws$ehc_surgeon_type1)),
  surgeon_type2  = summarize_posterior(exp(draws$ehc_surgeon_type2))
)

draws$contrast_ref_type1 <-
  draws$`b_log_ehc_z:social_type1`

draws$contrast_ref_type2 <-
  draws$`b_log_ehc_z:social_type2`

draws$contrast_surgeon_type1 <-
  draws$`b_log_ehc_z:social_type1` +
  draws$`b_log_ehc_z:subjectSurgeonfish:social_type1`

draws$contrast_surgeon_type2 <-
  draws$`b_log_ehc_z:social_type2` +
  draws$`b_log_ehc_z:subjectSurgeonfish:social_type2`

summarize_posterior(draws$contrast_ref_type1)
summarize_posterior(draws$contrast_ref_type2)
summarize_posterior(draws$contrast_surgeon_type1)
summarize_posterior(draws$contrast_surgeon_type2)

summarize_posterior(exp(draws$contrast_ref_type1))
summarize_posterior(exp(draws$contrast_ref_type2))
summarize_posterior(exp(draws$contrast_surgeon_type1))
summarize_posterior(exp(draws$contrast_surgeon_type2))


draws <- as_draws_df(bite_type_model)

draws$ehc_ref_type0 <-
  draws$b_log_ehc_z

draws$ehc_ref_type1 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:social_type1`

draws$ehc_ref_type2 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:social_type2`

draws$ehc_surgeon_type0 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:subjectSurgeonfish`

draws$ehc_surgeon_type1 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:subjectSurgeonfish` +
  draws$`b_log_ehc_z:social_type1` +
  draws$`b_log_ehc_z:subjectSurgeonfish:social_type1`

draws$ehc_surgeon_type2 <-
  draws$b_log_ehc_z +
  draws$`b_log_ehc_z:subjectSurgeonfish` +
  draws$`b_log_ehc_z:social_type2` +
  draws$`b_log_ehc_z:subjectSurgeonfish:social_type2`

rbind(
  ref_type0     = summarize_posterior(draws$ehc_ref_type0),
  ref_type1     = summarize_posterior(draws$ehc_ref_type1),
  ref_type2     = summarize_posterior(draws$ehc_ref_type2),
  surgeon_type0 = summarize_posterior(draws$ehc_surgeon_type0),
  surgeon_type1 = summarize_posterior(draws$ehc_surgeon_type1),
  surgeon_type2 = summarize_posterior(draws$ehc_surgeon_type2)
)
