# ========================================================================================= #
# PYI Impact Study                                                                          #
# Model outcome: in_housing_spell_on_date_18 using nn matched sample                        #
# Author: David Taylor                                                                      #
# Date: 20/01/2025                                                                          #
# ========================================================================================= #

# This code: 
#   Is divided into four sections;
#   1. Estimate Treatment Effects (ATT)
#   2. Subgroup analysis x male 
#   3. Subgroup analysis x aboriginal
#   4. Subgroup analysis x vulnerability   
#   5. Sensitivity analysis using tipping point analysis 
#   In each section, summary output is exported that allows for consolidation of results with other model output
#   
#   Notes:
#   * average potential outcomes and ATT and CATT estimates in RD, RR and OR are estimated using marginaleffects()
#   * a custom function is used to estimate cluster bootstrap standard errors and 95% confidence intervals

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages 
library(tidyverse)
library(MatchIt)
library(marginaleffects)
library(brglm2)
library(boot)
library(tipr)

# set up workspace preferences using custom helper function
workspace_setup()

# set seed
set.seed(268)

# read modelling data
modelling_data_file_location <- "P:/pyi/pyi_update/data/processed_data/pyi_modelling_data_nn_match.RDS" 
modelling_data <- readRDS(modelling_data_file_location)

# read sensitivity analysis inputs
tip_analysis_inputs_file_location <- "U:/pyi-paper/output/sensitivity_analysis/tip_analysis_inputs.RDS"
tip_analysis_inputs <- readRDS(tip_analysis_inputs_file_location)

#-------------------------------------------------------------------------------
# 1. Estimate Treatment Effects 
#-------------------------------------------------------------------------------

# in_housing_spell_on_date_18: ATT model
in_housing_spell_on_date_18_model <- glm(
  in_housing_spell_on_date_18 ~ 
    pyi_flag  + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data,
  family = quasibinomial(),
  weights = weights
)

# in_housing_spell_on_date_18: get average estimated potential outcomes in natural scale (i.e., probabilities)
in_housing_spell_on_date_18_model_potential_outcomes <- avg_predictions(
  in_housing_spell_on_date_18_model,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1)
)

# in_housing_spell_on_date_18: get ATT in RD
in_housing_spell_on_date_18_model_att_rd <- avg_comparisons(
  in_housing_spell_on_date_18_model,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1)
)

# in_housing_spell_on_date_18: get ATT in RR
in_housing_spell_on_date_18_model_att_rr <- avg_comparisons(
  in_housing_spell_on_date_18_model,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  comparison = "lnratioavg",
  transform = "exp"
)

# in_housing_spell_on_date_18: get ATT in OR
in_housing_spell_on_date_18_model_att_or <- avg_comparisons(
  in_housing_spell_on_date_18_model,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  comparison = "lnoravg",
  transform = "exp"
)

# bootstrap cluster se and confidence intervals for TE
#-----------------------------------------------------

# get unique pair ids from matched data
pair_ids = levels(modelling_data$subclass)

# create split indicies
split_inds <- split(
  seq_len(
    nrow(
      modelling_data)
  ),
  modelling_data$subclass
)

# declare bootstrap function for RD, RR and OR
cluster_boot_function <- function(pairs, i) {
  # extract units corresponding to selected pairs
  ids <- unlist(split_inds[pairs[i]])
  
  # subset model data with block bootstrapped indices
  boot_model_data <- modelling_data[ids, ]
  
  # refit model with bootstrapped data
  boot_model <- glm(
    in_housing_spell_on_date_18 ~ 
      pyi_flag  + 
      eligible_residential_care_placement +
      eligible_time_in_care +                                
      eligible_placement_instability +
      eligible_permanent_placement,
    data = boot_model_data, 
    family = quasibinomial(),
    weights = weights
  )
  
  # g-computation: subset treated units for att
  outcome_data <- subset(
    boot_model_data, 
    pyi_flag == 1)
  
  # estimate potential outcomes under intervention
  potential_outcomes_intervention <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data,
      pyi_flag = 1)
  )
  
  # estimate potential outcomes under comparison
  potential_outcomes_comparison <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data,
      pyi_flag = 0)
  )
  
  # estimate mean potential outcomes
  estimated_potential_outcomes_intervention <- mean(potential_outcomes_intervention)
  estimated_potential_outcomes_comparison <- mean(potential_outcomes_comparison) 
  
  # estimate logit-based smd
  log_or <- log(estimated_potential_outcomes_intervention/(1-estimated_potential_outcomes_intervention)) - log(estimated_potential_outcomes_comparison/(1-estimated_potential_outcomes_comparison))
  standardised_mean_difference <- (sqrt(3)/pi) * log_or
  
  # summary information
  risk_difference <- estimated_potential_outcomes_intervention - estimated_potential_outcomes_comparison
  relative_risk <- estimated_potential_outcomes_intervention / estimated_potential_outcomes_comparison
  odds_ratio <- (estimated_potential_outcomes_intervention / (1-estimated_potential_outcomes_intervention)) / 
    (estimated_potential_outcomes_comparison / (1-estimated_potential_outcomes_comparison))
  
  # return information
  return(c(
    rd = risk_difference,
    rr = relative_risk,
    or = odds_ratio,
    smd = standardised_mean_difference)
  )
}

# run bootstrap for all effects
in_housing_spell_on_date_18_att_estimates_boot <- boot(
  pair_ids,
  cluster_boot_function,
  R = 4999
)

# get smd
in_housing_spell_on_date_18_att_smd <- in_housing_spell_on_date_18_att_estimates_boot$t0 |>
  as_tibble() |>
  filter(
    row_number() == 4
  ) |>
  as.numeric()

# get boot se
in_housing_spell_on_date_18_att_bootstrap_se <- in_housing_spell_on_date_18_att_estimates_boot$t |>
  as_tibble() |>
  rename(
    "rd" = 1, 
    "rr" = 2,  
    "or" = 3,
    "smd" = 4
  ) |>
  summarise(
    across(
      everything(),
      \(x) sd(x),
      .names = "se_{.col}"
    )
  )

# get confidence intervals
in_housing_spell_on_date_18_model_att_ci_boot_rd <- boot.ci(
  in_housing_spell_on_date_18_att_estimates_boot, 
  type = "bca", 
  index = 1)

in_housing_spell_on_date_18_model_att_ci_boot_rr <- boot.ci(
  in_housing_spell_on_date_18_att_estimates_boot, 
  type = "bca", 
  index = 2)

in_housing_spell_on_date_18_model_att_ci_boot_or <- boot.ci(
  in_housing_spell_on_date_18_att_estimates_boot, 
  type = "bca", 
  index = 3)

in_housing_spell_on_date_18_model_att_ci_boot_smd <- boot.ci(
  in_housing_spell_on_date_18_att_estimates_boot, 
  type = "bca", 
  index = 4)

# extract information for export
in_housing_spell_on_date_18_model_family <- "Quasibinomial"
in_housing_spell_on_date_18_model_link_function <- "Logit"
in_housing_spell_on_date_18_model_num_obvs <- nobs(in_housing_spell_on_date_18_model)
in_housing_spell_on_date_18_model_potential_outcomes_intervention <- round(in_housing_spell_on_date_18_model_potential_outcomes$estimate[2], 3)
in_housing_spell_on_date_18_model_potential_outcomes_comparison <- round(in_housing_spell_on_date_18_model_potential_outcomes$estimate[1], 3)
in_housing_spell_on_date_18_model_att_rd_estimate <- round(in_housing_spell_on_date_18_model_att_rd$estimate, 3)
in_housing_spell_on_date_18_model_att_rd_se <- paste("(", round(in_housing_spell_on_date_18_att_bootstrap_se$se_rd, 3), ")", sep = "")
in_housing_spell_on_date_18_model_att_rd_conf_int <- paste("[", round(in_housing_spell_on_date_18_model_att_ci_boot_rd$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_att_ci_boot_rd$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_att_rr_estimate <- round(in_housing_spell_on_date_18_model_att_rr$estimate, 3)
in_housing_spell_on_date_18_model_att_rr_se <- paste("(", round(in_housing_spell_on_date_18_att_bootstrap_se$se_rr, 3), ")", sep = "")
in_housing_spell_on_date_18_model_att_rr_conf_int <- paste("[", round(in_housing_spell_on_date_18_model_att_rr$conf.low, 3), ", ", round(in_housing_spell_on_date_18_model_att_rr$conf.high, 3), "]", sep = "")
in_housing_spell_on_date_18_model_att_or_estimate <- round(in_housing_spell_on_date_18_model_att_or$estimate, 3)
in_housing_spell_on_date_18_model_att_or_se <- paste("(", round(in_housing_spell_on_date_18_att_bootstrap_se$se_or, 3), ")", sep = "")
in_housing_spell_on_date_18_model_att_or_conf_int <- paste("[", round(in_housing_spell_on_date_18_model_att_ci_boot_or$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_att_ci_boot_or$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_att_smd_estimate <- round(in_housing_spell_on_date_18_att_smd, 3)
in_housing_spell_on_date_18_model_att_smd_se <- paste("(", round(in_housing_spell_on_date_18_att_bootstrap_se$se_smd, 3), ")", sep = "")
in_housing_spell_on_date_18_model_att_smd_conf_int <- paste("[", round(in_housing_spell_on_date_18_model_att_ci_boot_smd$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_att_ci_boot_smd$bca[5], 3), "]", sep = "")

# bring together results in a summary table
in_housing_spell_on_date_18_results <- tibble(
  reported_metrics = c(
    "Model family",
    "Link function",
    "Estimated Potential Outcomes:",
    "Intervention group",
    "Comparison group",
    "Treatment Effects:",
    "Risk Difference:",
    "ATT estimate",
    "Standard error",
    "95% CI",
    "Relative Risk:",
    "ATT estimate",
    "Standard error",
    "95% CI",
    "Odds Ratio:",
    "ATT estimate",
    "Standard error",
    "95% CI",
    "Mean Difference:",
    "ATT estimate",
    "Standard error",
    "95% CI",
    "Standardised Mean Difference:",
    "ATT estimate",
    "Standard error",
    "95% CI",
    "Number of Observations"
  ),
  in_housing_spell_on_date_18 = c(
    in_housing_spell_on_date_18_model_family,
    in_housing_spell_on_date_18_model_link_function,
    NA,
    in_housing_spell_on_date_18_model_potential_outcomes_intervention,
    in_housing_spell_on_date_18_model_potential_outcomes_comparison,
    NA,
    NA,
    in_housing_spell_on_date_18_model_att_rd_estimate,
    in_housing_spell_on_date_18_model_att_rd_se,
    in_housing_spell_on_date_18_model_att_rd_conf_int,
    NA,
    in_housing_spell_on_date_18_model_att_rr_estimate,
    in_housing_spell_on_date_18_model_att_rr_se,
    in_housing_spell_on_date_18_model_att_rr_conf_int,
    NA,
    in_housing_spell_on_date_18_model_att_or_estimate,
    in_housing_spell_on_date_18_model_att_or_se,
    in_housing_spell_on_date_18_model_att_or_conf_int,
    NA,
    NA,
    NA,
    NA,
    NA,
    in_housing_spell_on_date_18_model_att_smd_estimate,
    in_housing_spell_on_date_18_model_att_smd_se,
    in_housing_spell_on_date_18_model_att_smd_conf_int,
    in_housing_spell_on_date_18_model_num_obvs)
)

# export treatment effect results
saveRDS(
  in_housing_spell_on_date_18_results,
  "U:/pyi-paper/output/treatment_effect_results/in_housing_spell_on_date_18_te_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 2. Subgroup analysis x Sex
#-------------------------------------------------------------------------------

# in_housing_spell_on_date_18: CATT model by sex
in_housing_spell_on_date_18_model_subgroup_sex <- glm(
  in_housing_spell_on_date_18 ~ 
    pyi_flag * 
    male,
  data = modelling_data,
  family = binomial(),
  weights = weights,
  method = "brglmFit",
  type = "MPL_Jeffreys"
)

# in_housing_spell_on_date_18: get average estimated potential outcomes by pyi_flag x sex 
in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes <- avg_predictions(
  in_housing_spell_on_date_18_model_subgroup_sex,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = c("pyi_flag", "male")
)

# in_housing_spell_on_date_18: get CATT by sex in RD
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_sex,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "male"
)

# in_housing_spell_on_date_18: get CATT by sex in RR
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_sex,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "male",
  comparison = "lnratioavg",
  transform = "exp"
)

# in_housing_spell_on_date_18: get CATT by sex in OR
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_sex,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "male",
  comparison = "lnoravg",
  transform = "exp"
)

# in_housing_spell_on_date_18: test for moderation of outcome by sex
in_housing_spell_on_date_18_model_subgroup_sex_moderation_test <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_sex,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "male",
  hypothesis = "pairwise"
)


# bootstrap cluster se and confidence intervals for subgroup analysis x sex
#----------------------------------------------------------------------------

# declare bootstrap function for RD, RR and OR
cluster_boot_function_subgroup_sex <- function(pairs, i) {
  # extract units corresponding to selected pairs
  ids <- unlist(split_inds[pairs[i]])
  
  # subset model data with block bootstrapped indices
  boot_model_data <- modelling_data[ids, ]
  
  # refit model with bootstrapped data
  boot_model <- glm(
    in_housing_spell_on_date_18 ~ 
      pyi_flag  * 
      male, 
    data = boot_model_data,
    family = binomial(),
    weights = weights,
    method = "brglmFit",
    type = "MPL_Jeffreys"
  )
  
  # g-computation: subset treated units for att
  
  # male subgroup
  outcome_data_male <- boot_model_data |>
    filter( 
      pyi_flag == 1,
      male == 1)
  
  # estimate potential outcomes under intervention conditional on being male
  potential_outcomes_intervention_male <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_male,
      pyi_flag = 1,
      male = 1)
  )
  
  # estimate potential outcomes under comparison conditional on being male
  potential_outcomes_comparison_male <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_male,
      pyi_flag = 0,
      male = 1
    )
  )
  
  # female subgroup
  outcome_data_female <- boot_model_data |>
    filter( 
      pyi_flag == 1,
      male == 0)
  
  # estimate potential outcomes under intervention conditional on being female
  potential_outcomes_intervention_female <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_female,
      pyi_flag = 1,
      male = 0)
  )
  
  # estimate potential outcomes under comparison conditional on being female
  potential_outcomes_comparison_female <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_female,
      pyi_flag = 0,
      male = 0)
  )
  
  # estimate mean potential outcomes
  estimated_potential_outcomes_intervention_male <- mean(potential_outcomes_intervention_male)
  estimated_potential_outcomes_intervention_female <- mean(potential_outcomes_intervention_female)
  estimated_potential_outcomes_comparison_male <- mean(potential_outcomes_comparison_male) 
  estimated_potential_outcomes_comparison_female <- mean(potential_outcomes_comparison_female) 
  
  # estimate logit-based smd
  log_or_male <- log(estimated_potential_outcomes_intervention_male/(1-estimated_potential_outcomes_intervention_male)) - log(estimated_potential_outcomes_comparison_male/(1-estimated_potential_outcomes_comparison_male))
  standardised_mean_difference_male <- (sqrt(3)/pi) * log_or_male
  
  log_or_female <- log(estimated_potential_outcomes_intervention_female/(1-estimated_potential_outcomes_intervention_female)) - log(estimated_potential_outcomes_comparison_female/(1-estimated_potential_outcomes_comparison_female))
  standardised_mean_difference_female <- (sqrt(3)/pi) * log_or_female
  
  # summary information
  risk_difference_male <- estimated_potential_outcomes_intervention_male - estimated_potential_outcomes_comparison_male
  risk_difference_female <- estimated_potential_outcomes_intervention_female - estimated_potential_outcomes_comparison_female
  relative_risk_male <- estimated_potential_outcomes_intervention_male / estimated_potential_outcomes_comparison_male
  relative_risk_female <- estimated_potential_outcomes_intervention_female / estimated_potential_outcomes_comparison_female
  odds_ratio_male <- (estimated_potential_outcomes_intervention_male / (1-estimated_potential_outcomes_intervention_male)) / 
    (estimated_potential_outcomes_comparison_male / (1-estimated_potential_outcomes_comparison_male))
  odds_ratio_female <- (estimated_potential_outcomes_intervention_female / (1-estimated_potential_outcomes_intervention_female)) / 
    (estimated_potential_outcomes_comparison_female / (1-estimated_potential_outcomes_comparison_female))
  
  # return information
  return(c(
    rd_male = risk_difference_male,
    rd_female = risk_difference_female,
    rr_male = relative_risk_male,
    rr_female = relative_risk_female,
    or_male = odds_ratio_male,
    or_female = odds_ratio_female,
    smd_male = standardised_mean_difference_male, 
    smd_female = standardised_mean_difference_female 
  ))
}

# run bootstrap for all effects
in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot <- boot(
  pair_ids,
  cluster_boot_function_subgroup_sex,
  R = 4999
)

# get smd
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_male <- in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot$t0 |>
  as_tibble() |>
  filter(
    row_number() == 7
  ) |>
  as.numeric()

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_female <- in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot$t0 |>
  as_tibble() |>
  filter(
    row_number() == 8
  ) |>
  as.numeric()

# get boot se
in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se <- in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot$t |>
  as_tibble() |>
  rename(
    "rd_male" = 1,
    "rd_female" = 2,
    "rr_male" = 3,
    "rr_female" = 4,
    "or_male" = 5,
    "or_female" = 6,
    "smd_male" = 7,
    "smd_female" = 8
  ) |>
  summarise(
    across(
      everything(),
      \(x) sd(x),
      .names = "se_{.col}"
    )
  )

# get confidence intervals
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rd_male <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 1)

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rd_female <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 2)

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rr_male <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 3)

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rr_female <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 4)

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_or_male <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 5)

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_or_female <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 6)

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_male <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 7)

in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_female <- boot.ci(
  in_housing_spell_on_date_18_subgroup_sex_conditional_att_estimates_boot, 
  type = "bca", 
  index = 8)

# extract information for export
in_housing_spell_on_date_18_model_subgroup_sex_family <- "Binomial (MPL Jeffreys bias-reduced)"
in_housing_spell_on_date_18_model_subgroup_sex_link_function <- "Logit"
in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_intervention_male <- round(in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes$estimate[4], 3)
in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_intervention_female <- round(in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes$estimate[3], 3)
in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_comparison_male <- round(in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_comparison_female <- round(in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_estimate_male <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_se_male <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_rd_male, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_conf_int_male <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rd_male$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rd_male$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_estimate_female <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_se_female <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_rd_female, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_conf_int_female <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rd_female$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rd_female$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_estimate_male <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_se_male <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_rr_male, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_conf_int_male <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rr_male$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rr_male$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_estimate_female <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_se_female <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_rr_female, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_conf_int_female <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rr_female$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_rr_female$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_estimate_male <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_se_male <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_or_male, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_conf_int_male <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_or_male$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_or_male$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_estimate_female <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_se_female <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_or_female, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_conf_int_female <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_or_female$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_or_female$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_male <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_male, 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_se_male <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_smd_male, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_conf_int_male <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_male$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_male$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_female <- round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_female, 3)
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_se_female <- paste("(", round(in_housing_spell_on_date_18_conditional_sex_att_bootstrap_se$se_smd_female, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_conf_int_female <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_female$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_female$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect <- round(in_housing_spell_on_date_18_model_subgroup_sex_moderation_test$estimate, 3)
in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_se <- paste("(", round(in_housing_spell_on_date_18_model_subgroup_sex_moderation_test$std.error, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_conf_int <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_sex_moderation_test$conf.low, 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_sex_moderation_test$conf.high, 3), "]", sep = "") 
in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_p_value <- round(in_housing_spell_on_date_18_model_subgroup_sex_moderation_test$p.value, 3)
in_housing_spell_on_date_18_model_subgroup_sex_num_obvs_male <- modelling_data |> filter(male == 1) |> nrow()
in_housing_spell_on_date_18_model_subgroup_sex_num_obvs_female <- modelling_data |> filter(male == 0) |> nrow()

# combine export information in data frame
in_housing_spell_on_date_18_subgroup_sex_results <- tibble(
  reported_metrics = c(
    "Model family",
    "Link function",
    "Estimated Potential Outcomes:",
    "Intervention group",
    "Comparison group",
    "Treatment Effects:",
    "Risk Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Relative Risk:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Odds Ratio:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Mean Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Standardised Mean Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Effect moderation:",
    "Difference",
    "Standard error",
    "95% CI",
    "p-value",
    "Number of Observations"
  ),
  in_housing_spell_on_date_18_male = c(
    in_housing_spell_on_date_18_model_subgroup_sex_family,
    in_housing_spell_on_date_18_model_subgroup_sex_link_function,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_intervention_male,
    in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_comparison_male,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_estimate_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_se_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_conf_int_male,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_estimate_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_se_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_conf_int_male,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_estimate_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_se_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_conf_int_male,
    NA,
    NA,
    NA,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_se_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_conf_int_male,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_se,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_conf_int,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_p_value,
    in_housing_spell_on_date_18_model_subgroup_sex_num_obvs_male
  ),
  in_housing_spell_on_date_18_female = c(
    in_housing_spell_on_date_18_model_subgroup_sex_family,
    in_housing_spell_on_date_18_model_subgroup_sex_link_function,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_intervention_female,
    in_housing_spell_on_date_18_model_subgroup_sex_potential_outcomes_comparison_female,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_estimate_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_se_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rd_conf_int_female,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_estimate_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_se_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_rr_conf_int_female,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_estimate_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_se_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_or_conf_int_female,
    NA,
    NA,
    NA,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_se_female,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_conf_int_female,
    NA,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_se,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_conf_int,
    in_housing_spell_on_date_18_model_subgroup_sex_moderated_effect_p_value,
    in_housing_spell_on_date_18_model_subgroup_sex_num_obvs_female
  ))

# export treatment effect results
saveRDS(
  in_housing_spell_on_date_18_subgroup_sex_results,
  "U:/pyi-paper/output/treatment_effect_results/in_housing_spell_on_date_18_subgroup_sex_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 2. Subgroup analysis x Aboriginal Status
#-------------------------------------------------------------------------------

# in_housing_spell_on_date_18: CATT model by aboriginal status
in_housing_spell_on_date_18_model_subgroup_aboriginal <- glm(
  in_housing_spell_on_date_18 ~ 
    pyi_flag  * 
    aboriginal,
  data = modelling_data,
  family = binomial(),
  weights = weights,
  method = "brglmFit",
  type = "MPL_Jeffreys"
)


# in_housing_spell_on_date_18: get average estimated potential outcomes by pyi_flag x aboriginal status 
in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes <- avg_predictions(
  in_housing_spell_on_date_18_model_subgroup_aboriginal,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = c("pyi_flag", "aboriginal")
)

# in_housing_spell_on_date_18: get CATT by aboriginal status in RD
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_aboriginal,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "aboriginal"
)

# in_housing_spell_on_date_18: get CATT by aboriginal status in RR
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_aboriginal,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "aboriginal",
  comparison = "lnratioavg",
  transform = "exp"
)

# in_housing_spell_on_date_18: get CATT by aboriginal status in OR
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_aboriginal,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "aboriginal",
  comparison = "lnoravg",
  transform = "exp"
)

# in_housing_spell_on_date_18: test for moderation of outcome by aboriginal status
in_housing_spell_on_date_18_model_subgroup_aboriginal_moderation_test <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_aboriginal,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "aboriginal",
  hypothesis = "pairwise"
)

# bootstrap cluster se and confidence intervals for subgroup analysis x Aboriginal status
#----------------------------------------------------------------------------------------

# declare bootstrap function for RD, RR and OR
cluster_boot_function_subgroup_aboriginal <- function(pairs, i) {
  # extract units corresponding to selected pairs
  ids <- unlist(split_inds[pairs[i]])
  
  # subset model data with block bootstrapped indices
  boot_model_data <- modelling_data[ids, ]
  
  # refit model with bootstrapped data
  # note: events are sparse in small sample so bias-reduced logistic regression is used
  boot_model <- glm(
    in_housing_spell_on_date_18 ~ 
      pyi_flag  * 
      aboriginal,
    data = boot_model_data,
    family = binomial(),
    weights = weights,
    method = "brglmFit",
    type = "MPL_Jeffreys")
  
  # g-computation: subset treated units for att
  
  # aboriginal subgroup
  outcome_data_aboriginal <- boot_model_data |>
    filter( 
      pyi_flag == 1,
      aboriginal == 1)
  
  # estimate potential outcomes under intervention conditional on being aboriginal
  potential_outcomes_intervention_aboriginal <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_aboriginal,
      pyi_flag = 1,
      aboriginal = 1)
  )
  
  # estimate potential outcomes under comparison conditional on being aboriginal
  potential_outcomes_comparison_aboriginal <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_aboriginal,
      pyi_flag = 0,
      aboriginal = 1)
  )
  
  # non_aboriginal subgroup
  outcome_data_non_aboriginal <- boot_model_data |>
    filter( 
      pyi_flag == 1,
      aboriginal == 0)
  
  # estimate potential outcomes under intervention conditional on being non_aboriginal
  potential_outcomes_intervention_non_aboriginal <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_non_aboriginal,
      pyi_flag = 1,
      aboriginal = 0)
  )
  
  # estimate potential outcomes under comparison conditional on being non_aboriginal
  potential_outcomes_comparison_non_aboriginal <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_non_aboriginal,
      pyi_flag = 0,
      aboriginal = 0)
  )
  
  # estimate mean potential outcomes
  estimated_potential_outcomes_intervention_aboriginal <- mean(potential_outcomes_intervention_aboriginal)
  estimated_potential_outcomes_intervention_non_aboriginal <- mean(potential_outcomes_intervention_non_aboriginal)
  estimated_potential_outcomes_comparison_aboriginal <- mean(potential_outcomes_comparison_aboriginal) 
  estimated_potential_outcomes_comparison_non_aboriginal <- mean(potential_outcomes_comparison_non_aboriginal) 
  
  # estimate logit-based smd
  log_or_aboriginal <- log(estimated_potential_outcomes_intervention_aboriginal/(1-estimated_potential_outcomes_intervention_aboriginal)) - log(estimated_potential_outcomes_comparison_aboriginal/(1-estimated_potential_outcomes_comparison_aboriginal))
  standardised_mean_difference_aboriginal <- (sqrt(3)/pi) * log_or_aboriginal
  
  log_or_non_aboriginal <- log(estimated_potential_outcomes_intervention_non_aboriginal/(1-estimated_potential_outcomes_intervention_non_aboriginal)) - log(estimated_potential_outcomes_comparison_non_aboriginal/(1-estimated_potential_outcomes_comparison_non_aboriginal))
  standardised_mean_difference_non_aboriginal <- (sqrt(3)/pi) * log_or_non_aboriginal
  
  # summary information
  risk_difference_aboriginal <- estimated_potential_outcomes_intervention_aboriginal - estimated_potential_outcomes_comparison_aboriginal
  risk_difference_non_aboriginal <- estimated_potential_outcomes_intervention_non_aboriginal - estimated_potential_outcomes_comparison_non_aboriginal
  relative_risk_aboriginal <- estimated_potential_outcomes_intervention_aboriginal / estimated_potential_outcomes_comparison_aboriginal
  relative_risk_non_aboriginal <- estimated_potential_outcomes_intervention_non_aboriginal / estimated_potential_outcomes_comparison_non_aboriginal
  odds_ratio_aboriginal <- (estimated_potential_outcomes_intervention_aboriginal / (1-estimated_potential_outcomes_intervention_aboriginal)) / 
    (estimated_potential_outcomes_comparison_aboriginal / (1-estimated_potential_outcomes_comparison_aboriginal))
  odds_ratio_non_aboriginal <- (estimated_potential_outcomes_intervention_non_aboriginal / (1-estimated_potential_outcomes_intervention_non_aboriginal)) / 
    (estimated_potential_outcomes_comparison_non_aboriginal / (1-estimated_potential_outcomes_comparison_non_aboriginal))
  
  # return information
  return(c(
    rd_aboriginal = risk_difference_aboriginal,
    rd_non_aboriginal = risk_difference_non_aboriginal,
    rr_aboriginal = relative_risk_aboriginal,
    rr_non_aboriginal = relative_risk_non_aboriginal,
    or_aboriginal = odds_ratio_aboriginal,
    or_non_aboriginal = odds_ratio_non_aboriginal,
    smd_aboriginal = standardised_mean_difference_aboriginal,
    smd_non_aboriginal = standardised_mean_difference_non_aboriginal)
  )
}

# run bootstrap for all effects
in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot <- boot(
  pair_ids,
  cluster_boot_function_subgroup_aboriginal,
  R = 4999
)

# get smd
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_aboriginal <- in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot$t0 |>
  as_tibble() |>
  filter(
    row_number() == 7
  ) |>
  as.numeric()

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_non_aboriginal <- in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot$t0 |>
  as_tibble() |>
  filter(
    row_number() == 8
  ) |>
  as.numeric()

# get boot se
in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se <- in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot$t |>
  as_tibble() |>
  rename(
    "rd_aboriginal" = 1,
    "rd_non_aboriginal" = 2,
    "rr_aboriginal" = 3,
    "rr_non_aboriginal" = 4,
    "or_aboriginal" = 5,
    "or_non_aboriginal" = 6,
    "smd_aboriginal" = 7,
    "smd_non_aboriginal" = 8
  ) |>
  summarise(
    across(
      everything(),
      \(x) sd(x),
      .names = "se_{.col}"
    )
  )

# get confidence intervals
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rd_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 1)

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rd_non_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 2)

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rr_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 3)

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rr_non_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 4)

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_or_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 5)

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_or_non_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 6)

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 7)

in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_non_aboriginal <- boot.ci(
  in_housing_spell_on_date_18_subgroup_aboriginal_conditional_att_estimates_boot, 
  type = "bca", 
  index = 8)

# extract information for export
in_housing_spell_on_date_18_model_subgroup_aboriginal_family <- "Binomial (MPL Jeffreys bias-reduced)"
in_housing_spell_on_date_18_model_subgroup_aboriginal_link_function <- "Logit"
in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_intervention_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes$estimate[4], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_intervention_non_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes$estimate[3], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_comparison_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_comparison_non_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_estimate_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_se_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_rd_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_conf_int_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rd_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rd_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_estimate_non_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_se_non_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_rd_non_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_conf_int_non_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rd_non_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rd_non_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_estimate_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_se_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_rr_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_conf_int_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rr_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rr_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_estimate_non_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_se_non_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_rr_non_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_conf_int_non_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rr_non_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_rr_non_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_estimate_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_se_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_or_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_conf_int_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_or_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_or_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_estimate_non_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_se_non_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_or_non_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_conf_int_non_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_or_non_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_or_non_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_aboriginal, 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_se_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_smd_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_conf_int_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_non_aboriginal <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_non_aboriginal, 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_se_non_aboriginal <- paste("(", round(in_housing_spell_on_date_18_conditional_aboriginal_att_bootstrap_se$se_smd_non_aboriginal, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_conf_int_non_aboriginal <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_non_aboriginal$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_non_aboriginal$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_moderation_test$estimate, 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_se <- paste("(", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_moderation_test$std.error, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_conf_int <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_moderation_test$conf.low, 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_aboriginal_moderation_test$conf.high, 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_p_value <- round(in_housing_spell_on_date_18_model_subgroup_aboriginal_moderation_test$p.value, 3)
in_housing_spell_on_date_18_model_subgroup_aboriginal_num_obvs_aboriginal <- modelling_data |> filter(aboriginal == 1) |> nrow()
in_housing_spell_on_date_18_model_subgroup_aboriginal_num_obvs_non_aboriginal <- modelling_data |> filter(aboriginal == 0) |> nrow()

# combine export information in data frame
in_housing_spell_on_date_18_subgroup_aboriginal_results <- tibble(
  reported_metrics = c(
    "Model family",
    "Link function",
    "Estimated Potential Outcomes:",
    "Intervention group",
    "Comparison group",
    "Treatment Effects:",
    "Risk Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Relative Risk:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Odds Ratio:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Mean Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Standardised Mean Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Effect moderation:",
    "Difference",
    "Standard error",
    "95% CI",
    "p-value",
    "Number of Observations"
  ),
  in_housing_spell_on_date_18_aboriginal = c(
    in_housing_spell_on_date_18_model_subgroup_aboriginal_family,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_link_function,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_intervention_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_comparison_aboriginal,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_estimate_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_se_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_conf_int_aboriginal,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_estimate_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_se_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_conf_int_aboriginal,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_estimate_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_se_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_conf_int_aboriginal,
    NA,
    NA,
    NA,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_se_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_conf_int_aboriginal,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_se,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_conf_int,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_p_value,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_num_obvs_aboriginal
  ),
  in_housing_spell_on_date_18_non_aboriginal = c(
    in_housing_spell_on_date_18_model_subgroup_aboriginal_family,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_link_function,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_intervention_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_potential_outcomes_comparison_non_aboriginal,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_estimate_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_se_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rd_conf_int_non_aboriginal,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_estimate_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_se_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_rr_conf_int_non_aboriginal,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_estimate_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_se_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_or_conf_int_non_aboriginal,
    NA,
    NA,
    NA,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_se_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_conf_int_non_aboriginal,
    NA,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_se,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_conf_int,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_moderated_effect_p_value,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_num_obvs_non_aboriginal
  ))

# export treatment effect results
saveRDS(
  in_housing_spell_on_date_18_subgroup_aboriginal_results,
  "U:/pyi-paper/output/treatment_effect_results/in_housing_spell_on_date_18_subgroup_aboriginal_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 4. Subgroup analysis x housing vulnerability
#-------------------------------------------------------------------------------

# in_housing_spell_on_date_18: CATT model by housing vulnerability
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability <- glm(
  in_housing_spell_on_date_18 ~ 
    pyi_flag * 
    any_housing_spell_before_18,
  data = modelling_data,
  family = binomial(),
  weights = weights,
  method = "brglmFit",
  type = "MPL_Jeffreys"
)

# in_housing_spell_on_date_18: get average estimated potential outcomes by pyi_flag x housing vulnerability 
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes <- avg_predictions(
  in_housing_spell_on_date_18_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = c("pyi_flag", "any_housing_spell_before_18")
)

# in_housing_spell_on_date_18: get CATT by housing vulnerability in RD
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "any_housing_spell_before_18"
)

# in_housing_spell_on_date_18: get CATT by housing vulnerability in RR
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "any_housing_spell_before_18",
  comparison = "lnratioavg",
  transform = "exp"
)

# in_housing_spell_on_date_18: get CATT by housing vulnerability in OR
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "any_housing_spell_before_18",
  comparison = "lnoravg",
  transform = "exp"
)

# in_housing_spell_on_date_18: test for moderation of outcome by housing vulnerability
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderation_test <- avg_comparisons(
  in_housing_spell_on_date_18_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "any_housing_spell_before_18",
  hypothesis = "pairwise"
)


# bootstrap cluster se and confidence intervals for subgroup analysis x housing vulnerability
#----------------------------------------------------------------------------

# declare bootstrap function for RD, RR and OR
cluster_boot_function_subgroup_housing_vulnerability <- function(pairs, i) {
  # extract units corresponding to selected pairs
  ids <- unlist(split_inds[pairs[i]])
  
  # subset model data with block bootstrapped indices
  boot_model_data <- modelling_data[ids, ]
  
  # refit model with bootstrapped data
  boot_model <- glm(
    in_housing_spell_on_date_18 ~ 
      pyi_flag  * 
      any_housing_spell_before_18, 
    data = boot_model_data,
    family = binomial(),
    weights = weights,
    method = "brglmFit",
    type = "MPL_Jeffreys"
  )
  
  # g-computation: subset treated units for att
  
  # prior_homelessness subgroup
  outcome_data_prior_homelessness <- boot_model_data |>
    filter( 
      pyi_flag == 1,
      any_housing_spell_before_18 == 1)
  
  # estimate potential outcomes under intervention conditional on prior homelessness
  potential_outcomes_intervention_prior_homelessness <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_prior_homelessness,
      pyi_flag = 1,
      any_housing_spell_before_18 = 1)
  )
  
  # estimate potential outcomes under comparison conditional on prior homelessness prior_homelessness
  potential_outcomes_comparison_prior_homelessness <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_prior_homelessness,
      pyi_flag = 0,
      any_housing_spell_before_18 = 1
    )
  )
  
  # no_prior_homelessness subgroup
  outcome_data_no_prior_homelessness <- boot_model_data |>
    filter( 
      pyi_flag == 1,
      any_housing_spell_before_18 == 0)
  
  # estimate potential outcomes under intervention conditional on no prior homelessness
  potential_outcomes_intervention_no_prior_homelessness <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_no_prior_homelessness,
      pyi_flag = 1,
      any_housing_spell_before_18 = 0)
  )
  
  # estimate potential outcomes under comparison conditional on no prior homelessness
  potential_outcomes_comparison_no_prior_homelessness <- predict(
    boot_model,
    type = "response",
    newdata = transform(
      outcome_data_no_prior_homelessness,
      pyi_flag = 0,
      any_housing_spell_before_18 = 0)
  )
  
  # estimate mean potential outcomes
  estimated_potential_outcomes_intervention_prior_homelessness <- mean(potential_outcomes_intervention_prior_homelessness)
  estimated_potential_outcomes_intervention_no_prior_homelessness <- mean(potential_outcomes_intervention_no_prior_homelessness)
  estimated_potential_outcomes_comparison_prior_homelessness <- mean(potential_outcomes_comparison_prior_homelessness) 
  estimated_potential_outcomes_comparison_no_prior_homelessness <- mean(potential_outcomes_comparison_no_prior_homelessness) 
  
  # estimate logit-based smd
  log_or_prior_homelessness <- log(estimated_potential_outcomes_intervention_prior_homelessness/(1-estimated_potential_outcomes_intervention_prior_homelessness)) - log(estimated_potential_outcomes_comparison_prior_homelessness/(1-estimated_potential_outcomes_comparison_prior_homelessness))
  standardised_mean_difference_prior_homelessness <- (sqrt(3)/pi) * log_or_prior_homelessness
  
  log_or_no_prior_homelessness <- log(estimated_potential_outcomes_intervention_no_prior_homelessness/(1-estimated_potential_outcomes_intervention_no_prior_homelessness)) - log(estimated_potential_outcomes_comparison_no_prior_homelessness/(1-estimated_potential_outcomes_comparison_no_prior_homelessness))
  standardised_mean_difference_no_prior_homelessness <- (sqrt(3)/pi) * log_or_no_prior_homelessness
  
  # summary information
  risk_difference_prior_homelessness <- estimated_potential_outcomes_intervention_prior_homelessness - estimated_potential_outcomes_comparison_prior_homelessness
  risk_difference_no_prior_homelessness <- estimated_potential_outcomes_intervention_no_prior_homelessness - estimated_potential_outcomes_comparison_no_prior_homelessness
  relative_risk_prior_homelessness <- estimated_potential_outcomes_intervention_prior_homelessness / estimated_potential_outcomes_comparison_prior_homelessness
  relative_risk_no_prior_homelessness <- estimated_potential_outcomes_intervention_no_prior_homelessness / estimated_potential_outcomes_comparison_no_prior_homelessness
  odds_ratio_prior_homelessness <- (estimated_potential_outcomes_intervention_prior_homelessness / (1-estimated_potential_outcomes_intervention_prior_homelessness)) / 
    (estimated_potential_outcomes_comparison_prior_homelessness / (1-estimated_potential_outcomes_comparison_prior_homelessness))
  odds_ratio_no_prior_homelessness <- (estimated_potential_outcomes_intervention_no_prior_homelessness / (1-estimated_potential_outcomes_intervention_no_prior_homelessness)) / 
    (estimated_potential_outcomes_comparison_no_prior_homelessness / (1-estimated_potential_outcomes_comparison_no_prior_homelessness))
  
  # return information
  return(c(
    rd_prior_homelessness = risk_difference_prior_homelessness,
    rd_no_prior_homelessness = risk_difference_no_prior_homelessness,
    rr_prior_homelessness = relative_risk_prior_homelessness,
    rr_no_prior_homelessness = relative_risk_no_prior_homelessness,
    or_prior_homelessness = odds_ratio_prior_homelessness,
    or_no_prior_homelessness = odds_ratio_no_prior_homelessness,
    smd_prior_homelessness = standardised_mean_difference_prior_homelessness, 
    smd_no_prior_homelessness = standardised_mean_difference_no_prior_homelessness 
  ))
}

# run bootstrap for all effects
in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot <- boot(
  pair_ids,
  cluster_boot_function_subgroup_housing_vulnerability,
  R = 4999
)

# get smd
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_prior_homelessness <- in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot$t0 |>
  as_tibble() |>
  filter(
    row_number() == 7
  ) |>
  as.numeric()

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_no_prior_homelessness <- in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot$t0 |>
  as_tibble() |>
  filter(
    row_number() == 8
  ) |>
  as.numeric()

# get boot se
in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se <- in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot$t |>
  as_tibble() |>
  rename(
    "rd_prior_homelessness" = 1,
    "rd_no_prior_homelessness" = 2,
    "rr_prior_homelessness" = 3,
    "rr_no_prior_homelessness" = 4,
    "or_prior_homelessness" = 5,
    "or_no_prior_homelessness" = 6,
    "smd_prior_homelessness" = 7,
    "smd_no_prior_homelessness" = 8
  ) |>
  summarise(
    across(
      everything(),
      \(x) sd(x),
      .names = "se_{.col}"
    )
  )

# get confidence intervals
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rd_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 1)

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rd_no_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 2)

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rr_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 3)

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rr_no_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 4)

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_or_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 5)

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_or_no_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 6)

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 7)

in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_no_prior_homelessness <- boot.ci(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_conditional_att_estimates_boot, 
  type = "bca", 
  index = 8)

# extract information for export
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_family <- "Binomial (MPL Jeffreys bias-reduced)"
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_link_function <- "Logit"
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_intervention_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes$estimate[4], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_intervention_no_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes$estimate[3], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_comparison_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_comparison_no_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_estimate_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_se_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_rd_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_conf_int_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rd_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rd_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_estimate_no_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_se_no_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_rd_no_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_conf_int_no_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rd_no_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rd_no_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_estimate_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_se_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_rr_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_conf_int_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rr_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rr_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_estimate_no_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_se_no_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_rr_no_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_conf_int_no_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rr_no_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_rr_no_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_estimate_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or$estimate[2], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_se_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_or_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_conf_int_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_or_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_or_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_estimate_no_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or$estimate[1], 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_se_no_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_or_no_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_conf_int_no_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_or_no_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_or_no_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_prior_homelessness, 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_se_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_smd_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_conf_int_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_no_prior_homelessness <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_no_prior_homelessness, 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_se_no_prior_homelessness <- paste("(", round(in_housing_spell_on_date_18_conditional_housing_vulnerability_att_bootstrap_se$se_smd_no_prior_homelessness, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_conf_int_no_prior_homelessness <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_no_prior_homelessness$bca[4], 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_no_prior_homelessness$bca[5], 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderation_test$estimate, 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_se <- paste("(", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderation_test$std.error, 3), ")", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_conf_int <- paste("[", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderation_test$conf.low, 3), ", ", round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderation_test$conf.high, 3), "]", sep = "")
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_p_value <- round(in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderation_test$p.value, 3)
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_num_obvs_prior_homelessness <- modelling_data |> filter(any_housing_spell_before_18 == 1) |> nrow()
in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_num_obvs_no_prior_homelessness <- modelling_data |> filter(any_housing_spell_before_18 == 0) |> nrow()

# combine export information in data frame
in_housing_spell_on_date_18_subgroup_housing_vulnerability_results <- tibble(
  reported_metrics = c(
    "Model family",
    "Link function",
    "Estimated Potential Outcomes:",
    "Intervention group",
    "Comparison group",
    "Treatment Effects:",
    "Risk Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Relative Risk:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Odds Ratio:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Mean Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Standardised Mean Difference:",
    "CATT estimate",
    "Standard error",
    "95% CI",
    "Effect moderation:",
    "Difference",
    "Standard error",
    "95% CI",
    "p-value",
    "Number of Observations"
  ),
  in_housing_spell_on_date_18_prior_homelessness = c(
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_family,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_link_function,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_intervention_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_comparison_prior_homelessness,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_estimate_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_se_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_conf_int_prior_homelessness,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_estimate_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_se_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_conf_int_prior_homelessness,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_estimate_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_se_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_conf_int_prior_homelessness,
    NA,
    NA,
    NA,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_se_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_conf_int_prior_homelessness,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_se,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_conf_int,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_p_value,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_num_obvs_prior_homelessness
  ),
  in_housing_spell_on_date_18_no_prior_homelessness = c(
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_family,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_link_function,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_intervention_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_potential_outcomes_comparison_no_prior_homelessness,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_estimate_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_se_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rd_conf_int_no_prior_homelessness,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_estimate_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_se_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_rr_conf_int_no_prior_homelessness,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_estimate_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_se_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_or_conf_int_no_prior_homelessness,
    NA,
    NA,
    NA,
    NA,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_se_no_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_conf_int_no_prior_homelessness,
    NA,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_se,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_conf_int,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_moderated_effect_p_value,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_num_obvs_no_prior_homelessness
  ))

# export treatment effect results
saveRDS(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_results,
  "U:/pyi-paper/output/treatment_effect_results/in_housing_spell_on_date_18_subgroup_housing_vulnerability_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 5. Sensitivity Analysis using Tipping Point Analysis 
#-------------------------------------------------------------------------------

# sensitivity analysis using tipping point analysis for ATT
in_housing_spell_on_date_18_tip_high_prevalence_scenario <- tipr::tip_with_binary(
  effect_observed = in_housing_spell_on_date_18_model_att_rr$estimate,
  exposed_confounder_prev = tip_analysis_inputs$prevalence_treated[1],
  unexposed_confounder_prev = tip_analysis_inputs$prevalence_untreated[1], 
  verbose = TRUE
)

in_housing_spell_on_date_18_tip_medium_prevalence_scenario <- tipr::tip_with_binary(
  effect_observed = in_housing_spell_on_date_18_model_att_rr$estimate,
  exposed_confounder_prev = tip_analysis_inputs$prevalence_treated[2],
  unexposed_confounder_prev = tip_analysis_inputs$prevalence_untreated[2], 
  verbose = TRUE
)

# note: this call is commented out because the following error term, stops the code running: 
# Given these prevalences (`unexposed_confounder_prev`: 0.05, `exposed_confounder_prev`: 0.25), there does not exist an unmeasured confounder that could tip this. 
#in_housing_spell_on_date_18_tip_low_prevalence_scenario <- tipr::tip_with_binary(
#  effect_observed = in_housing_spell_on_date_18_model_att_rr$estimate,
#  exposed_confounder_prev = tip_analysis_inputs$prevalence_treated[3],
#  unexposed_confounder_prev = tip_analysis_inputs$prevalence_untreated[3], 
#  verbose = TRUE
#) 

# combine and export sensitivity analysis results
in_housing_spell_on_date_18_tip_results <- tibble(
  scenario = c(
    "High", 
    "Medium", 
    "Low"),
  confounder_prevalence_exposed = tip_analysis_inputs$prevalence_treated,
  confounder_prevalence_unexposed = tip_analysis_inputs$prevalence_untreated,
  confounder_outcome_relationship_in_housing_spell_on_date_18 = c(
    in_housing_spell_on_date_18_tip_high_prevalence_scenario$confounder_outcome_effect,
    in_housing_spell_on_date_18_tip_medium_prevalence_scenario$confounder_outcome_effect,
    NA)
)

# export sensitivity analysis results
saveRDS(
  in_housing_spell_on_date_18_tip_results,
  "U:/pyi-paper/output/sensitivity_analysis/in_housing_spell_on_date_18_tip_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 6. Export SMD values for plotting 
#-------------------------------------------------------------------------------

in_housing_spell_on_date_18_smd_results <- tibble(
  outcome = c(
    "In homelessness spell on 18th birthday",
    "In homelessness spell on 18th birthday",
    "In homelessness spell on 18th birthday",
    "In homelessness spell on 18th birthday",
    "In homelessness spell on 18th birthday",
    "In homelessness spell on 18th birthday",
    "In homelessness spell on 18th birthday"
  ),
  group = c(
    "Overall",
    "Sex",
    "Sex",
    "Aboriginal",
    "Aboriginal",
    "Housing vulnerability",
    "Housing vulnerability"
  ),
  strata = c(
    "Overall",
    "Male",
    "Female",
    "Aboriginal",
    "Non-Aboriginal",
    "Prior homelessness",
    "No prior homelessness"
  ),
  smd = c(
    in_housing_spell_on_date_18_att_smd,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_male,
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_smd_estimate_female,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_smd_estimate_non_aboriginal,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_prior_homelessness,
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_smd_estimate_no_prior_homelessness
  ),
  ci_low = c(
    in_housing_spell_on_date_18_model_att_ci_boot_smd$bca[4],
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_male$bca[4],
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_female$bca[4],    
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_aboriginal$bca[4],
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_non_aboriginal$bca[4],
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_prior_homelessness$bca[4],
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_no_prior_homelessness$bca[4]
  ),
  ci_high = c(
    in_housing_spell_on_date_18_model_att_ci_boot_smd$bca[5],
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_male$bca[5],
    in_housing_spell_on_date_18_model_subgroup_sex_conditional_att_ci_boot_smd_female$bca[5],
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_aboriginal$bca[5],
    in_housing_spell_on_date_18_model_subgroup_aboriginal_conditional_att_ci_boot_smd_non_aboriginal$bca[5],
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_prior_homelessness$bca[5],
    in_housing_spell_on_date_18_model_subgroup_housing_vulnerability_conditional_att_ci_boot_smd_no_prior_homelessness$bca[5]
  ))


# export smd results
saveRDS(
  in_housing_spell_on_date_18_smd_results,
  "U:/pyi-paper/output/plot_data/in_housing_spell_on_date_18_smd_results_nn_match.RDS"
)
