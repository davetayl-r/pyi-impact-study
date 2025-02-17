# ========================================================================================= #
# PYI Impact Study                                                                          #
# Model outcome: count_total_spells_18_19 using nn matched sample                           #
# Author: David Taylor                                                                      #
# Date: 23/01/2025                                                                          #
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
library(MASS)
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

# A. Poisson model
#-------------------------------------------------------------------------------

# count_total_spells_18_19: ATT model
count_total_spells_18_19_poisson_model <- glm(
  count_total_spells_18_19 ~ 
    pyi_flag  + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data,
  weights = weights,
  family = poisson(link = "log")
)

# count_total_spells_18_19: get ATT
count_total_spells_18_19_poisson_model_att <- avg_comparisons(
  count_total_spells_18_19_poisson_model,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1)
)

# B. Negative binomial model
#-------------------------------------------------------------------------------
count_total_spells_18_19_neg_binom_model <- MASS::glm.nb(
  count_total_spells_18_19 ~ 
    pyi_flag  + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data,
  weights = weights
)

count_total_spells_18_19_neg_binom_model_att <- avg_comparisons(
  count_total_spells_18_19_neg_binom_model,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1)
)

## NOTE: 
# Results from the poisson and negative binomial model are functionally the same, prefer the poisson model 

# count_total_spells_18_19: get average estimated potential outcomes in natural scale (i.e., means)
count_total_spells_18_19_poisson_model_potential_outcomes <- avg_predictions(
  count_total_spells_18_19_poisson_model,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1)
)

# get sample size in each group
count_total_spells_18_19_intervention_n <- modelling_data |>
  filter(
    pyi_flag == 1
  ) |>
  nrow()

count_total_spells_18_19_comparison_n <- modelling_data |>
  filter(
    pyi_flag == 0
  ) |>
  nrow()

# estimate sd in each group
count_total_spells_18_19_intervention_sd <- modelling_data |>
  filter(
    pyi_flag == 1
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_comparison_sd <- modelling_data |>
  filter(
    pyi_flag == 0
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

# estimate pooled sd
count_total_spells_18_19_pooled_sd <- sqrt(
  (count_total_spells_18_19_intervention_sd^2 + count_total_spells_18_19_comparison_sd^2) / 2
)

# estimate smd
count_total_spells_18_19_poisson_model_att_smd <- (count_total_spells_18_19_poisson_model_att$estimate / count_total_spells_18_19_pooled_sd) |>
  as.numeric()

# estimate smd se
count_total_spells_18_19_poisson_model_att_smd_se <- sqrt(
  (count_total_spells_18_19_intervention_n + count_total_spells_18_19_comparison_n) / (count_total_spells_18_19_intervention_n * count_total_spells_18_19_comparison_n) +
    count_total_spells_18_19_poisson_model_att_smd^2/(2*(count_total_spells_18_19_intervention_n + count_total_spells_18_19_comparison_n))  
)

# estimate smd confidence interval
count_total_spells_18_19_poisson_model_att_smd_ci_low <- count_total_spells_18_19_poisson_model_att_smd - 1.96 * count_total_spells_18_19_poisson_model_att_smd_se
count_total_spells_18_19_poisson_model_att_smd_ci_high <- count_total_spells_18_19_poisson_model_att_smd + 1.96 * count_total_spells_18_19_poisson_model_att_smd_se

# extract information for export
count_total_spells_18_19_poisson_model_family <- "Poisson"
count_total_spells_18_19_poisson_model_link_function <- "Log"
count_total_spells_18_19_num_obvs <- nobs(count_total_spells_18_19_poisson_model)
count_total_spells_18_19_potential_outcomes_intervention <- round(count_total_spells_18_19_poisson_model_potential_outcomes$estimate[2], 3)
count_total_spells_18_19_potential_outcomes_comparison <- round(count_total_spells_18_19_poisson_model_potential_outcomes$estimate[1], 3)
count_total_spells_18_19_att_mean_estimate <- round(count_total_spells_18_19_poisson_model_att$estimate, 3)
count_total_spells_18_19_att_mean_se <- paste("(", round(count_total_spells_18_19_poisson_model_att$std.error, 3), ")", sep = "")
count_total_spells_18_19_att_mean_conf_int <- paste("[", round(count_total_spells_18_19_poisson_model_att$conf.low, 3), ", ", round(count_total_spells_18_19_poisson_model_att$conf.high, 3), "]", sep = "")
count_total_spells_18_19_att_smd_estimate <- round(count_total_spells_18_19_poisson_model_att_smd, 3)
count_total_spells_18_19_att_smd_se <- paste("(", round(count_total_spells_18_19_poisson_model_att_smd_se, 3), ")", sep = "")
count_total_spells_18_19_att_smd_conf_int <- paste("[", round(count_total_spells_18_19_poisson_model_att_smd_ci_low, 3), ", ", round(count_total_spells_18_19_poisson_model_att_smd_ci_high, 3), "]", sep = "")


count_total_spells_18_19_results <- tibble(
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
  count_total_spells_18_19 = c(
    count_total_spells_18_19_poisson_model_family,
    count_total_spells_18_19_poisson_model_link_function,
    NA,
    count_total_spells_18_19_potential_outcomes_intervention,
    count_total_spells_18_19_potential_outcomes_comparison,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    count_total_spells_18_19_att_mean_estimate,
    count_total_spells_18_19_att_mean_se,
    count_total_spells_18_19_att_mean_conf_int,
    NA,
    count_total_spells_18_19_att_smd_estimate,
    count_total_spells_18_19_att_smd_se,
    count_total_spells_18_19_att_smd_conf_int,
    count_total_spells_18_19_num_obvs)
)

# export treatment effect results
saveRDS(
  count_total_spells_18_19_results,
  "U:/pyi-paper/output/treatment_effect_results/count_total_spells_18_19_te_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 2. Subgroup analysis x Sex
#-------------------------------------------------------------------------------

# A. Poisson model
#-------------------------------------------------------------------------------

# count_total_spells_18_19: CATT model by sex
count_total_spells_18_19_poisson_model_subgroup_sex <- glm(
  count_total_spells_18_19 ~ 
    pyi_flag  * 
    male,
  data = modelling_data,
  weights = weights,
  family = poisson(link = "log")
)

# count_total_spells_18_19: get CATT in means
count_total_spells_18_19_poisson_model_subgroup_sex_means <- avg_comparisons(
  count_total_spells_18_19_poisson_model_subgroup_sex,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "male"
)

# B. Negative binomial model
#-------------------------------------------------------------------------------

# count_total_spells_18_19: CATT model by sex
count_total_spells_18_19_neg_binom_model_subgroup_sex <- MASS::glm.nb(
  count_total_spells_18_19 ~ 
    pyi_flag  *
    male,
  data = modelling_data,
  weights = weights
)

# count_total_spells_18_19: get CATT in means
count_total_spells_18_19_neg_binom_model_subgroup_sex_means <- avg_comparisons(
  count_total_spells_18_19_neg_binom_model_subgroup_sex,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "male"
)

## NOTE: 
# Results from the poisson and negative binomial model are functionally the same, prefer the poisson model 

# count_total_spells_18_19: get average estimated potential outcomes by pyi_flag x sex 
count_total_spells_18_19_poisson_model_subgroup_sex_potential_outcomes <- avg_predictions(
  count_total_spells_18_19_poisson_model_subgroup_sex,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = c("pyi_flag", "male")
)

# count_total_spells_18_19: test for moderation of outcome by sex
count_total_spells_18_19_poisson_model_subgroup_sex_moderation_test <- avg_comparisons(
  count_total_spells_18_19_poisson_model_subgroup_sex,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "male",
  hypothesis = "pairwise"
)

# get sample size in each subgroup
count_total_spells_18_19_intervention_male_n <- modelling_data |>
  filter(
    pyi_flag == 1,
    male == 1
  ) |>
  nrow()

count_total_spells_18_19_intervention_female_n <- modelling_data |>
  filter(
    pyi_flag == 1,
    male == 0
  ) |>
  nrow()

count_total_spells_18_19_comparison_male_n <- modelling_data |>
  filter(
    pyi_flag == 0,
    male == 1
  ) |>
  nrow()

count_total_spells_18_19_comparison_female_n <- modelling_data |>
  filter(
    pyi_flag == 0,
    male == 0
  ) |>
  nrow()

# estimate sd in each subgroup
count_total_spells_18_19_intervention_male_sd <- modelling_data |>
  filter(
    pyi_flag == 1,
    male == 1
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_comparison_male_sd <- modelling_data |>
  filter(
    pyi_flag == 0,
    male == 1
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_intervention_female_sd <- modelling_data |>
  filter(
    pyi_flag == 1,
    male == 0
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_comparison_female_sd <- modelling_data |>
  filter(
    pyi_flag == 0,
    male == 0
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

# estimate pooled sd for each subgroup
count_total_spells_18_19_male_pooled_sd <- sqrt(
  (count_total_spells_18_19_intervention_male_sd^2 + count_total_spells_18_19_comparison_male_sd^2) / 2
)

count_total_spells_18_19_female_pooled_sd <- sqrt(
  (count_total_spells_18_19_intervention_female_sd^2 + count_total_spells_18_19_comparison_female_sd^2) / 2
)

# estimate smd for each subgroup
count_total_spells_18_19_male_smd <- (count_total_spells_18_19_poisson_model_subgroup_sex_means$estimate[2] / count_total_spells_18_19_male_pooled_sd) |>
  as.numeric()

count_total_spells_18_19_female_smd <- (count_total_spells_18_19_poisson_model_subgroup_sex_means$estimate[1] / count_total_spells_18_19_female_pooled_sd) |>
  as.numeric()


# estimate smd se for each subgroup
count_total_spells_18_19_model_male_smd_se <- sqrt(
  (count_total_spells_18_19_intervention_male_n + count_total_spells_18_19_comparison_male_n) / (count_total_spells_18_19_intervention_male_n * count_total_spells_18_19_comparison_male_n) +
    count_total_spells_18_19_male_smd^2/(2*(count_total_spells_18_19_intervention_male_n + count_total_spells_18_19_comparison_male_n))  
)

count_total_spells_18_19_model_female_smd_se <- sqrt(
  (count_total_spells_18_19_intervention_female_n + count_total_spells_18_19_comparison_female_n) / (count_total_spells_18_19_intervention_female_n * count_total_spells_18_19_comparison_female_n) +
    count_total_spells_18_19_female_smd^2/(2*(count_total_spells_18_19_intervention_female_n + count_total_spells_18_19_comparison_female_n))  
)

# estimate smd confidence interval for each subgroup
count_total_spells_18_19_model_male_smd_ci_low <- count_total_spells_18_19_male_smd - 1.96 * count_total_spells_18_19_model_male_smd_se
count_total_spells_18_19_model_male_smd_ci_high <- count_total_spells_18_19_male_smd + 1.96 * count_total_spells_18_19_model_male_smd_se
count_total_spells_18_19_model_female_smd_ci_low <- count_total_spells_18_19_female_smd - 1.96 * count_total_spells_18_19_model_female_smd_se
count_total_spells_18_19_model_female_smd_ci_high <- count_total_spells_18_19_female_smd + 1.96 * count_total_spells_18_19_model_female_smd_se

# extract information for export
count_total_spells_18_19_model_family <- "Poisson"
count_total_spells_18_19_model_link_function <- "Log"
count_total_spells_18_19_model_subgroup_sex_potential_outcomes_intervention_male <- round(count_total_spells_18_19_poisson_model_subgroup_sex_potential_outcomes$estimate[4], 3)
count_total_spells_18_19_model_subgroup_sex_potential_outcomes_intervention_female <- round(count_total_spells_18_19_poisson_model_subgroup_sex_potential_outcomes$estimate[3], 3)
count_total_spells_18_19_model_subgroup_sex_potential_outcomes_comparison_male <- round(count_total_spells_18_19_poisson_model_subgroup_sex_potential_outcomes$estimate[2], 3)
count_total_spells_18_19_model_subgroup_sex_potential_outcomes_comparison_female <- round(count_total_spells_18_19_poisson_model_subgroup_sex_potential_outcomes$estimate[1], 3)
count_total_spells_18_19_model_poisson_subgroup_sex_estimate_male <- round(count_total_spells_18_19_poisson_model_subgroup_sex_means$estimate[2], 3)
count_total_spells_18_19_mean_estimate_male <- round(count_total_spells_18_19_poisson_model_subgroup_sex_means$estimate[2], 3)
count_total_spells_18_19_mean_se_male <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_sex_means$std.error[2], 3), ")", sep = "")
count_total_spells_18_19_mean_conf_int_male <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_sex_means$conf.low[2], 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_sex_means$conf.high[2], 3), "]", sep = "")
count_total_spells_18_19_mean_estimate_female <- round(count_total_spells_18_19_poisson_model_subgroup_sex_means$estimate[1], 3)
count_total_spells_18_19_mean_se_female <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_sex_means$std.error[1], 3), ")", sep = "")
count_total_spells_18_19_mean_conf_int_female <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_sex_means$conf.low[1], 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_sex_means$conf.high[1], 3), "]", sep = "")
count_total_spells_18_19_smd_estimate_male <- round(count_total_spells_18_19_male_smd, 3)
count_total_spells_18_19_smd_se_male <- paste("(", round(count_total_spells_18_19_model_male_smd_se, 3), ")", sep = "")
count_total_spells_18_19_smd_conf_int_male <- paste("[", round(count_total_spells_18_19_model_male_smd_ci_low, 3), ", ", round(count_total_spells_18_19_model_male_smd_ci_high, 3), "]", sep = "")
count_total_spells_18_19_smd_estimate_female <- round(count_total_spells_18_19_female_smd, 3)
count_total_spells_18_19_smd_se_female <- paste("(", round(count_total_spells_18_19_model_female_smd_se, 3), ")", sep = "")
count_total_spells_18_19_smd_conf_int_female <- paste("[", round(count_total_spells_18_19_model_female_smd_ci_low, 3), ", ", round(count_total_spells_18_19_model_female_smd_ci_high, 3), "]", sep = "")
count_total_spells_18_19_model_subgroup_sex_moderated_effect <- round(count_total_spells_18_19_poisson_model_subgroup_sex_moderation_test$estimate, 3)
count_total_spells_18_19_model_subgroup_sex_moderated_effect_se <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_sex_moderation_test$std.error, 3), ")", sep = "")
count_total_spells_18_19_model_subgroup_sex_moderated_effect_conf_int <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_sex_moderation_test$conf.low, 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_sex_moderation_test$conf.high, 3), "]", sep = "")
count_total_spells_18_19_model_subgroup_sex_moderated_effect_p_value <- round(count_total_spells_18_19_poisson_model_subgroup_sex_moderation_test$p.value, 3)
count_total_spells_18_19_model_subgroup_sex_num_obvs_male <- modelling_data |> filter(male == 1) |> nrow()
count_total_spells_18_19_model_subgroup_sex_num_obvs_female <- modelling_data |> filter(male == 0) |> nrow()

count_total_spells_18_19_subgroup_sex_results <- tibble(
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
  count_total_spells_18_19_male = c(
    count_total_spells_18_19_model_family,
    count_total_spells_18_19_model_link_function,
    NA,
    count_total_spells_18_19_model_subgroup_sex_potential_outcomes_intervention_male,
    count_total_spells_18_19_model_subgroup_sex_potential_outcomes_comparison_male,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    count_total_spells_18_19_mean_estimate_male,
    count_total_spells_18_19_mean_se_male,
    count_total_spells_18_19_mean_conf_int_male,
    NA,
    count_total_spells_18_19_smd_estimate_male,
    count_total_spells_18_19_smd_se_male,
    count_total_spells_18_19_smd_conf_int_male,
    NA,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect_se,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect_conf_int,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect_p_value,
    count_total_spells_18_19_model_subgroup_sex_num_obvs_male
  ),
  count_total_spells_18_19_female = c(
    count_total_spells_18_19_model_family,
    count_total_spells_18_19_model_link_function,
    NA,
    count_total_spells_18_19_model_subgroup_sex_potential_outcomes_intervention_female,
    count_total_spells_18_19_model_subgroup_sex_potential_outcomes_comparison_female,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    count_total_spells_18_19_mean_estimate_female,
    count_total_spells_18_19_mean_se_female,
    count_total_spells_18_19_mean_conf_int_female,
    NA,
    count_total_spells_18_19_smd_estimate_female,
    count_total_spells_18_19_smd_se_female,
    count_total_spells_18_19_smd_conf_int_female,
    NA,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect_se,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect_conf_int,
    count_total_spells_18_19_model_subgroup_sex_moderated_effect_p_value,
    count_total_spells_18_19_model_subgroup_sex_num_obvs_female)
)

# export treatment effect results
saveRDS(
  count_total_spells_18_19_subgroup_sex_results,
  "U:/pyi-paper/output/treatment_effect_results/count_total_spells_18_19_subgroup_sex_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 2. Subgroup analysis x Aboriginal Status
#-------------------------------------------------------------------------------

# A. Poisson model
#-------------------------------------------------------------------------------

# count_total_spells_18_19: CATT model by aboriginal status
count_total_spells_18_19_poisson_model_subgroup_aboriginal <- glm(
  count_total_spells_18_19 ~ 
    pyi_flag  * 
    aboriginal,
  data = modelling_data,
  weights = weights,
  family = poisson(link = "log")
)

# count_total_spells_18_19: get CATT in means
count_total_spells_18_19_poisson_model_subgroup_aboriginal_means <- avg_comparisons(
  count_total_spells_18_19_poisson_model_subgroup_aboriginal,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "aboriginal"
)

# B. Negative binomial model
#-------------------------------------------------------------------------------

# count_total_spells_18_19: CATT model by aboriginal status
count_total_spells_18_19_neg_binom_model_subgroup_aboriginal <- MASS::glm.nb(
  count_total_spells_18_19 ~ 
    pyi_flag  *
    aboriginal,
  data = modelling_data,
  weights = weights
)

# count_total_spells_18_19: get CATT in means
count_total_spells_18_19_neg_binom_model_subgroup_aboriginal_means <- avg_comparisons(
  count_total_spells_18_19_neg_binom_model_subgroup_aboriginal,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "aboriginal"
)

## NOTE: 
# Results from the poisson and negative binomial model are functionally the same, prefer the poisson model 

# count_total_spells_18_19: get average estimated potential outcomes by pyi_flag x aboriginal status 
count_total_spells_18_19_poisson_model_subgroup_aboriginal_potential_outcomes <- avg_predictions(
  count_total_spells_18_19_poisson_model_subgroup_aboriginal,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = c("pyi_flag", "aboriginal")
)

# count_total_spells_18_19: test for moderation of outcome by aboriginal status
count_total_spells_18_19_poisson_model_subgroup_aboriginal_moderation_test <- avg_comparisons(
  count_total_spells_18_19_poisson_model_subgroup_aboriginal,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "aboriginal",
  hypothesis = "pairwise"
)

# get sample size in each subgroup
count_total_spells_18_19_intervention_aboriginal_n <- modelling_data |>
  filter(
    pyi_flag == 1,
    aboriginal == 1
  ) |>
  nrow()

count_total_spells_18_19_intervention_non_aboriginal_n <- modelling_data |>
  filter(
    pyi_flag == 1,
    aboriginal == 0
  ) |>
  nrow()

count_total_spells_18_19_comparison_aboriginal_n <- modelling_data |>
  filter(
    pyi_flag == 0,
    aboriginal == 1
  ) |>
  nrow()

count_total_spells_18_19_comparison_non_aboriginal_n <- modelling_data |>
  filter(
    pyi_flag == 0,
    aboriginal == 0
  ) |>
  nrow()

# estimate sd in each subgroup
count_total_spells_18_19_intervention_aboriginal_sd <- modelling_data |>
  filter(
    pyi_flag == 1,
    aboriginal == 1
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_comparison_aboriginal_sd <- modelling_data |>
  filter(
    pyi_flag == 0,
    aboriginal == 1
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_intervention_non_aboriginal_sd <- modelling_data |>
  filter(
    pyi_flag == 1,
    aboriginal == 0
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_comparison_non_aboriginal_sd <- modelling_data |>
  filter(
    pyi_flag == 0,
    aboriginal == 0
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

# estimate pooled sd for each subgroup
count_total_spells_18_19_aboriginal_pooled_sd <- sqrt(
  (count_total_spells_18_19_intervention_aboriginal_sd^2 + count_total_spells_18_19_comparison_aboriginal_sd^2) / 2
)

count_total_spells_18_19_non_aboriginal_pooled_sd <- sqrt(
  (count_total_spells_18_19_intervention_non_aboriginal_sd^2 + count_total_spells_18_19_comparison_non_aboriginal_sd^2) / 2
)

# estimate smd for each subgroup
count_total_spells_18_19_aboriginal_smd <- (count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$estimate[2] / count_total_spells_18_19_aboriginal_pooled_sd) |>
  as.numeric()

count_total_spells_18_19_non_aboriginal_smd <- (count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$estimate[1] / count_total_spells_18_19_non_aboriginal_pooled_sd) |>
  as.numeric()


# estimate smd se for each subgroup
count_total_spells_18_19_model_aboriginal_smd_se <- sqrt(
  (count_total_spells_18_19_intervention_aboriginal_n + count_total_spells_18_19_comparison_aboriginal_n) / (count_total_spells_18_19_intervention_aboriginal_n * count_total_spells_18_19_comparison_aboriginal_n) +
    count_total_spells_18_19_aboriginal_smd^2/(2*(count_total_spells_18_19_intervention_aboriginal_n + count_total_spells_18_19_comparison_aboriginal_n))  
)

count_total_spells_18_19_model_non_aboriginal_smd_se <- sqrt(
  (count_total_spells_18_19_intervention_non_aboriginal_n + count_total_spells_18_19_comparison_non_aboriginal_n) / (count_total_spells_18_19_intervention_non_aboriginal_n * count_total_spells_18_19_comparison_non_aboriginal_n) +
    count_total_spells_18_19_non_aboriginal_smd^2/(2*(count_total_spells_18_19_intervention_non_aboriginal_n + count_total_spells_18_19_comparison_non_aboriginal_n))  
)

# estimate smd confidence interval for each subgroup
count_total_spells_18_19_model_aboriginal_smd_ci_low <- count_total_spells_18_19_aboriginal_smd - 1.96 * count_total_spells_18_19_model_aboriginal_smd_se
count_total_spells_18_19_model_aboriginal_smd_ci_high <- count_total_spells_18_19_aboriginal_smd + 1.96 * count_total_spells_18_19_model_aboriginal_smd_se
count_total_spells_18_19_model_non_aboriginal_smd_ci_low <- count_total_spells_18_19_non_aboriginal_smd - 1.96 * count_total_spells_18_19_model_non_aboriginal_smd_se
count_total_spells_18_19_model_non_aboriginal_smd_ci_high <- count_total_spells_18_19_non_aboriginal_smd + 1.96 * count_total_spells_18_19_model_non_aboriginal_smd_se

# extract information for export
count_total_spells_18_19_model_family <- "Poisson"
count_total_spells_18_19_model_link_function <- "Log"
count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_intervention_aboriginal <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_potential_outcomes$estimate[4], 3)
count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_intervention_non_aboriginal <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_potential_outcomes$estimate[3], 3)
count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_comparison_aboriginal <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_potential_outcomes$estimate[2], 3)
count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_comparison_non_aboriginal <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_potential_outcomes$estimate[1], 3)
count_total_spells_18_19_model_poisson_subgroup_aboriginal_estimate_aboriginal <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$estimate[2], 3)
count_total_spells_18_19_mean_estimate_aboriginal <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$estimate[2], 3)
count_total_spells_18_19_mean_se_aboriginal <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$std.error[2], 3), ")", sep = "")
count_total_spells_18_19_mean_conf_int_aboriginal <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$conf.low[2], 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$conf.high[2], 3), "]", sep = "")
count_total_spells_18_19_mean_estimate_non_aboriginal <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$estimate[1], 3)
count_total_spells_18_19_mean_se_non_aboriginal <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$std.error[1], 3), ")", sep = "")
count_total_spells_18_19_mean_conf_int_non_aboriginal <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$conf.low[1], 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_means$conf.high[1], 3), "]", sep = "")
count_total_spells_18_19_smd_estimate_aboriginal <- round(count_total_spells_18_19_aboriginal_smd, 3)
count_total_spells_18_19_smd_se_aboriginal <- paste("(", round(count_total_spells_18_19_model_aboriginal_smd_se, 3), ")", sep = "")
count_total_spells_18_19_smd_conf_int_aboriginal <- paste("[", round(count_total_spells_18_19_model_aboriginal_smd_ci_low, 3), ", ", round(count_total_spells_18_19_model_aboriginal_smd_ci_high, 3), "]", sep = "")
count_total_spells_18_19_smd_estimate_non_aboriginal <- round(count_total_spells_18_19_non_aboriginal_smd, 3)
count_total_spells_18_19_smd_se_non_aboriginal <- paste("(", round(count_total_spells_18_19_model_non_aboriginal_smd_se, 3), ")", sep = "")
count_total_spells_18_19_smd_conf_int_non_aboriginal <- paste("[", round(count_total_spells_18_19_model_non_aboriginal_smd_ci_low, 3), ", ", round(count_total_spells_18_19_model_non_aboriginal_smd_ci_high, 3), "]", sep = "")
count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_moderation_test$estimate, 3)
count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_se <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_moderation_test$std.error, 3), ")", sep = "")
count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_conf_int <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_moderation_test$conf.low, 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_moderation_test$conf.high, 3), "]", sep = "")
count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_p_value <- round(count_total_spells_18_19_poisson_model_subgroup_aboriginal_moderation_test$p.value, 3)
count_total_spells_18_19_model_subgroup_aboriginal_num_obvs_aboriginal <- modelling_data |> filter(aboriginal == 1) |> nrow()
count_total_spells_18_19_model_subgroup_aboriginal_num_obvs_non_aboriginal <- modelling_data |> filter(aboriginal == 0) |> nrow()

count_total_spells_18_19_subgroup_aboriginal_results <- tibble(
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
  count_total_spells_18_19_aboriginal = c(
    count_total_spells_18_19_model_family,
    count_total_spells_18_19_model_link_function,
    NA,
    count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_intervention_aboriginal,
    count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_comparison_aboriginal,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    count_total_spells_18_19_mean_estimate_aboriginal,
    count_total_spells_18_19_mean_se_aboriginal,
    count_total_spells_18_19_mean_conf_int_aboriginal,
    NA,
    count_total_spells_18_19_smd_estimate_aboriginal,
    count_total_spells_18_19_smd_se_aboriginal,
    count_total_spells_18_19_smd_conf_int_aboriginal,
    NA,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_se,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_conf_int,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_p_value,
    count_total_spells_18_19_model_subgroup_aboriginal_num_obvs_aboriginal
  ),
  count_total_spells_18_19_non_aboriginal = c(
    count_total_spells_18_19_model_family,
    count_total_spells_18_19_model_link_function,
    NA,
    count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_intervention_non_aboriginal,
    count_total_spells_18_19_model_subgroup_aboriginal_potential_outcomes_comparison_non_aboriginal,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    count_total_spells_18_19_mean_estimate_non_aboriginal,
    count_total_spells_18_19_mean_se_non_aboriginal,
    count_total_spells_18_19_mean_conf_int_non_aboriginal,
    NA,
    count_total_spells_18_19_smd_estimate_non_aboriginal,
    count_total_spells_18_19_smd_se_non_aboriginal,
    count_total_spells_18_19_smd_conf_int_non_aboriginal,
    NA,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_se,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_conf_int,
    count_total_spells_18_19_model_subgroup_aboriginal_moderated_effect_p_value,
    count_total_spells_18_19_model_subgroup_aboriginal_num_obvs_non_aboriginal)
)

# export treatment effect results
saveRDS(
  count_total_spells_18_19_subgroup_aboriginal_results,
  "U:/pyi-paper/output/treatment_effect_results/count_total_spells_18_19_subgroup_aboriginal_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 4. Subgroup analysis x housing vulnerability
#-------------------------------------------------------------------------------

# A. Poisson model
#-------------------------------------------------------------------------------

# count_total_spells_18_19: CATT model by housing vulnerability
count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability <- glm(
  count_total_spells_18_19 ~ 
    pyi_flag  * 
    any_housing_spell_before_18,
  data = modelling_data,
  weights = weights,
  family = poisson(link = "log")
)

# count_total_spells_18_19: get CATT in means
count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means <- avg_comparisons(
  count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "any_housing_spell_before_18"
)

# B. Negative binomial model
#-------------------------------------------------------------------------------

# count_total_spells_18_19: CATT model by housing vulnerability
count_total_spells_18_19_neg_binom_model_subgroup_housing_vulnerability <- MASS::glm.nb(
  count_total_spells_18_19 ~ 
    pyi_flag  *
    any_housing_spell_before_18,
  data = modelling_data,
  weights = weights
)

# count_total_spells_18_19: get CATT in means
count_total_spells_18_19_neg_binom_model_subgroup_housing_vulnerability_means <- avg_comparisons(
  count_total_spells_18_19_neg_binom_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "any_housing_spell_before_18"
)

## NOTE: 
# Results from the poisson and negative binomial model are functionally the same, prefer the poisson model 

# count_total_spells_18_19: get average estimated potential outcomes by pyi_flag x housing vulnerability 
count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_potential_outcomes <- avg_predictions(
  count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  comparison = "noravg",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = c("pyi_flag", "any_housing_spell_before_18")
)

# count_total_spells_18_19: test for moderation of outcome by housing vulnerability
count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_moderation_test <- avg_comparisons(
  count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability,
  variables = "pyi_flag",
  vcov = ~subclass,
  newdata = subset(pyi_flag == 1),
  by = "any_housing_spell_before_18",
  hypothesis = "pairwise"
)

# get sample size in each subgroup
count_total_spells_18_19_intervention_prior_homelessness_n <- modelling_data |>
  filter(
    pyi_flag == 1,
    any_housing_spell_before_18 == 1
  ) |>
  nrow()

count_total_spells_18_19_intervention_no_prior_homelessness_n <- modelling_data |>
  filter(
    pyi_flag == 1,
    any_housing_spell_before_18 == 0
  ) |>
  nrow()

count_total_spells_18_19_comparison_prior_homelessness_n <- modelling_data |>
  filter(
    pyi_flag == 0,
    any_housing_spell_before_18 == 1
  ) |>
  nrow()

count_total_spells_18_19_comparison_no_prior_homelessness_n <- modelling_data |>
  filter(
    pyi_flag == 0,
    any_housing_spell_before_18 == 0
  ) |>
  nrow()

# estimate sd in each subgroup
count_total_spells_18_19_intervention_prior_homelessness_sd <- modelling_data |>
  filter(
    pyi_flag == 1,
    any_housing_spell_before_18 == 1
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_comparison_prior_homelessness_sd <- modelling_data |>
  filter(
    pyi_flag == 0,
    any_housing_spell_before_18 == 1
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_intervention_no_prior_homelessness_sd <- modelling_data |>
  filter(
    pyi_flag == 1,
    any_housing_spell_before_18 == 0
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

count_total_spells_18_19_comparison_no_prior_homelessness_sd <- modelling_data |>
  filter(
    pyi_flag == 0,
    any_housing_spell_before_18 == 0
  ) |>
  select(
    count_total_spells_18_19
  ) |>
  summarise(
    sd = sd(count_total_spells_18_19)
  )

# estimate pooled sd for each subgroup
count_total_spells_18_19_prior_homelessness_pooled_sd <- sqrt(
  (count_total_spells_18_19_intervention_prior_homelessness_sd^2 + count_total_spells_18_19_comparison_prior_homelessness_sd^2) / 2
)

count_total_spells_18_19_no_prior_homelessness_pooled_sd <- sqrt(
  (count_total_spells_18_19_intervention_no_prior_homelessness_sd^2 + count_total_spells_18_19_comparison_no_prior_homelessness_sd^2) / 2
)

# estimate smd for each subgroup
count_total_spells_18_19_prior_homelessness_smd <- (count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$estimate[2] / count_total_spells_18_19_prior_homelessness_pooled_sd) |>
  as.numeric()

count_total_spells_18_19_no_prior_homelessness_smd <- (count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$estimate[1] / count_total_spells_18_19_no_prior_homelessness_pooled_sd) |>
  as.numeric()


# estimate smd se for each subgroup
count_total_spells_18_19_model_prior_homelessness_smd_se <- sqrt(
  (count_total_spells_18_19_intervention_prior_homelessness_n + count_total_spells_18_19_comparison_prior_homelessness_n) / (count_total_spells_18_19_intervention_prior_homelessness_n * count_total_spells_18_19_comparison_prior_homelessness_n) +
    count_total_spells_18_19_prior_homelessness_smd^2/(2*(count_total_spells_18_19_intervention_prior_homelessness_n + count_total_spells_18_19_comparison_prior_homelessness_n))  
)

count_total_spells_18_19_model_no_prior_homelessness_smd_se <- sqrt(
  (count_total_spells_18_19_intervention_no_prior_homelessness_n + count_total_spells_18_19_comparison_no_prior_homelessness_n) / (count_total_spells_18_19_intervention_no_prior_homelessness_n * count_total_spells_18_19_comparison_no_prior_homelessness_n) +
    count_total_spells_18_19_no_prior_homelessness_smd^2/(2*(count_total_spells_18_19_intervention_no_prior_homelessness_n + count_total_spells_18_19_comparison_no_prior_homelessness_n))  
)

# estimate smd confidence interval for each subgroup
count_total_spells_18_19_model_prior_homelessness_smd_ci_low <- count_total_spells_18_19_prior_homelessness_smd - 1.96 * count_total_spells_18_19_model_prior_homelessness_smd_se
count_total_spells_18_19_model_prior_homelessness_smd_ci_high <- count_total_spells_18_19_prior_homelessness_smd + 1.96 * count_total_spells_18_19_model_prior_homelessness_smd_se
count_total_spells_18_19_model_no_prior_homelessness_smd_ci_low <- count_total_spells_18_19_no_prior_homelessness_smd - 1.96 * count_total_spells_18_19_model_no_prior_homelessness_smd_se
count_total_spells_18_19_model_no_prior_homelessness_smd_ci_high <- count_total_spells_18_19_no_prior_homelessness_smd + 1.96 * count_total_spells_18_19_model_no_prior_homelessness_smd_se

# extract information for export
count_total_spells_18_19_model_family <- "Poisson"
count_total_spells_18_19_model_link_function <- "Log"
count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_intervention_prior_homelessness <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_potential_outcomes$estimate[4], 3)
count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_intervention_no_prior_homelessness <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_potential_outcomes$estimate[3], 3)
count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_comparison_prior_homelessness <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_potential_outcomes$estimate[2], 3)
count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_comparison_no_prior_homelessness <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_potential_outcomes$estimate[1], 3)
count_total_spells_18_19_model_poisson_subgroup_housing_vulnerability_estimate_prior_homelessness <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$estimate[2], 3)
count_total_spells_18_19_mean_estimate_prior_homelessness <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$estimate[2], 3)
count_total_spells_18_19_mean_se_prior_homelessness <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$std.error[2], 3), ")", sep = "")
count_total_spells_18_19_mean_conf_int_prior_homelessness <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$conf.low[2], 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$conf.high[2], 3), "]", sep = "")
count_total_spells_18_19_mean_estimate_no_prior_homelessness <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$estimate[1], 3)
count_total_spells_18_19_mean_se_no_prior_homelessness <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$std.error[1], 3), ")", sep = "")
count_total_spells_18_19_mean_conf_int_no_prior_homelessness <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$conf.low[1], 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_means$conf.high[1], 3), "]", sep = "")
count_total_spells_18_19_smd_estimate_prior_homelessness <- round(count_total_spells_18_19_prior_homelessness_smd, 3)
count_total_spells_18_19_smd_se_prior_homelessness <- paste("(", round(count_total_spells_18_19_model_prior_homelessness_smd_se, 3), ")", sep = "")
count_total_spells_18_19_smd_conf_int_prior_homelessness <- paste("[", round(count_total_spells_18_19_model_prior_homelessness_smd_ci_low, 3), ", ", round(count_total_spells_18_19_model_prior_homelessness_smd_ci_high, 3), "]", sep = "")
count_total_spells_18_19_smd_estimate_no_prior_homelessness <- round(count_total_spells_18_19_no_prior_homelessness_smd, 3)
count_total_spells_18_19_smd_se_no_prior_homelessness <- paste("(", round(count_total_spells_18_19_model_no_prior_homelessness_smd_se, 3), ")", sep = "")
count_total_spells_18_19_smd_conf_int_no_prior_homelessness <- paste("[", round(count_total_spells_18_19_model_no_prior_homelessness_smd_ci_low, 3), ", ", round(count_total_spells_18_19_model_no_prior_homelessness_smd_ci_high, 3), "]", sep = "")
count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_moderation_test$estimate, 3)
count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_se <- paste("(", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_moderation_test$std.error, 3), ")", sep = "")
count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_conf_int <- paste("[", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_moderation_test$conf.low, 3), ", ", round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_moderation_test$conf.high, 3), "]", sep = "")
count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_p_value <- round(count_total_spells_18_19_poisson_model_subgroup_housing_vulnerability_moderation_test$p.value, 3)
count_total_spells_18_19_model_subgroup_housing_vulnerability_num_obvs_prior_homelessness <- modelling_data |> filter(any_housing_spell_before_18 == 1) |> nrow()
count_total_spells_18_19_model_subgroup_housing_vulnerability_num_obvs_no_prior_homelessness <- modelling_data |> filter(any_housing_spell_before_18 == 0) |> nrow()

count_total_spells_18_19_subgroup_housing_vulnerability_results <- tibble(
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
  count_total_spells_18_19_prior_homelessness = c(
    count_total_spells_18_19_model_family,
    count_total_spells_18_19_model_link_function,
    NA,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_intervention_prior_homelessness,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_comparison_prior_homelessness,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    count_total_spells_18_19_mean_estimate_prior_homelessness,
    count_total_spells_18_19_mean_se_prior_homelessness,
    count_total_spells_18_19_mean_conf_int_prior_homelessness,
    NA,
    count_total_spells_18_19_smd_estimate_prior_homelessness,
    count_total_spells_18_19_smd_se_prior_homelessness,
    count_total_spells_18_19_smd_conf_int_prior_homelessness,
    NA,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_se,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_conf_int,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_p_value,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_num_obvs_prior_homelessness
  ),
  count_total_spells_18_19_no_prior_homelessness = c(
    count_total_spells_18_19_model_family,
    count_total_spells_18_19_model_link_function,
    NA,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_intervention_no_prior_homelessness,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_potential_outcomes_comparison_no_prior_homelessness,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    NA,
    count_total_spells_18_19_mean_estimate_no_prior_homelessness,
    count_total_spells_18_19_mean_se_no_prior_homelessness,
    count_total_spells_18_19_mean_conf_int_no_prior_homelessness,
    NA,
    count_total_spells_18_19_smd_estimate_no_prior_homelessness,
    count_total_spells_18_19_smd_se_no_prior_homelessness,
    count_total_spells_18_19_smd_conf_int_no_prior_homelessness,
    NA,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_se,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_conf_int,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_moderated_effect_p_value,
    count_total_spells_18_19_model_subgroup_housing_vulnerability_num_obvs_no_prior_homelessness)
)

# export treatment effect results
saveRDS(
  count_total_spells_18_19_subgroup_housing_vulnerability_results,
  "U:/pyi-paper/output/treatment_effect_results/count_total_spells_18_19_subgroup_housing_vulnerability_results_nn_match.RDS"
)

#-------------------------------------------------------------------------------
# 5. Sensitivity Analysis using Tipping Point Analysis 
#-------------------------------------------------------------------------------

# not done for continuous outcomes

#-------------------------------------------------------------------------------
# 6. Export SMD values for plotting 
#-------------------------------------------------------------------------------

count_total_spells_18_19_smd_results <- tibble(
  outcome = c(
    "Number of distinct homelessness spells between 18th and 19th birthday",
    "Number of distinct homelessness spells between 18th and 19th birthday",
    "Number of distinct homelessness spells between 18th and 19th birthday",
    "Number of distinct homelessness spells between 18th and 19th birthday",
    "Number of distinct homelessness spells between 18th and 19th birthday",
    "Number of distinct homelessness spells between 18th and 19th birthday",
    "Number of distinct homelessness spells between 18th and 19th birthday"
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
    count_total_spells_18_19_att_smd_estimate,
    count_total_spells_18_19_male_smd,
    count_total_spells_18_19_female_smd,
    count_total_spells_18_19_aboriginal_smd,
    count_total_spells_18_19_non_aboriginal_smd,
    count_total_spells_18_19_prior_homelessness_smd,
    count_total_spells_18_19_no_prior_homelessness_smd
  ),
  ci_low = c(
    count_total_spells_18_19_poisson_model_att_smd_ci_low,
    count_total_spells_18_19_model_male_smd_ci_low,
    count_total_spells_18_19_model_female_smd_ci_low,    
    count_total_spells_18_19_model_aboriginal_smd_ci_low,
    count_total_spells_18_19_model_non_aboriginal_smd_ci_low,
    count_total_spells_18_19_model_prior_homelessness_smd_ci_low,
    count_total_spells_18_19_model_no_prior_homelessness_smd_ci_low
  ),
  ci_high = c(
    count_total_spells_18_19_poisson_model_att_smd_ci_high,
    count_total_spells_18_19_model_male_smd_ci_high,
    count_total_spells_18_19_model_female_smd_ci_high,
    count_total_spells_18_19_model_aboriginal_smd_ci_high,
    count_total_spells_18_19_model_non_aboriginal_smd_ci_high,
    count_total_spells_18_19_model_prior_homelessness_smd_ci_high,
    count_total_spells_18_19_model_no_prior_homelessness_smd_ci_high
  ))


# export smd results
saveRDS(
  count_total_spells_18_19_smd_results,
  "U:/pyi-paper/output/plot_data/count_total_spells_18_19_smd_results_nn_match.RDS"
)
