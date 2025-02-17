# ========================================================================= #
# PYI Impact Study                                                          #
# Plot TE Heterogeneity: in_housing_spell_on_date_19                        #
# Author: David Taylor                                                      #
# Date: 28/01/2025                                                          #
# ========================================================================= #

# This code examines treatment effect heterogeneity for each matching specification and exports the results for visualisation

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

# read modelling data
modelling_data_nn_match_file_location <- "P:/pyi/pyi_update/data/processed_data/pyi_modelling_data_nn_match.RDS" 
modelling_data_nn_match <- readRDS(modelling_data_nn_match_file_location)

modelling_data_full_match_file_location <- "P:/pyi/pyi_update/data/processed_data/pyi_modelling_data_full_match.RDS" 
modelling_data_full_match <- readRDS(modelling_data_full_match_file_location)

#-------------------------------------------------------------------------------
# 1. Nearest Neighbour Matching
#-------------------------------------------------------------------------------

# ATT model
in_housing_spell_on_date_19_model_att_nn_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag  + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_nn_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x sex
in_housing_spell_on_date_19_model_catt_sex_nn_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * male + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_nn_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x aboriginal status
in_housing_spell_on_date_19_model_catt_aboriginal_nn_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * aboriginal + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_nn_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x housing vulnerability
in_housing_spell_on_date_19_model_catt_housing_vulnerability_nn_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * any_housing_spell_before_18 + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_nn_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x residential care
in_housing_spell_on_date_19_model_catt_residential_care_nn_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * eligible_residential_care_placement + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_nn_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x placement instability
in_housing_spell_on_date_19_model_catt_placement_instability_nn_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * eligible_placement_instability + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_nn_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x kinship care
in_housing_spell_on_date_19_model_catt_kinship_care_nn_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * eligible_kinship_care_placement + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_nn_match,
  family = quasibinomial(),
  weights = weights
)

# Step B: Estimate Individual TE using G-computation
#---------------------------------------------------

treatment_effect_heterogeneity_in_housing_spell_on_date_19_nn_match <- modelling_data_nn_match |>
  mutate(
    # estimate potential outcomes under intervention for ATT
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention = predict(
      in_housing_spell_on_date_19_model_att_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for ATT
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison = predict(
      in_housing_spell_on_date_19_model_att_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 0)
    ),
    # estimate ATT
    in_housing_spell_on_date_19_model_nn_match_treatment_effect = in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention - in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison,
    # estimate potential outcomes under intervention for CATT x sex
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_sex = predict(
      in_housing_spell_on_date_19_model_catt_sex_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT sex
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_sex = predict(
      in_housing_spell_on_date_19_model_catt_sex_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 0)
    ),
    # estimate CATT sex
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_sex = in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_sex - in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_sex,
    # estimate potential outcomes under intervention for CATT x aboriginal status
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_aboriginal = predict(
      in_housing_spell_on_date_19_model_catt_sex_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x aboriginal status
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_aboriginal = predict(
      in_housing_spell_on_date_19_model_catt_sex_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 0)
    ),
    # estimate CATT sex
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_aboriginal = in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_aboriginal - in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_aboriginal,
    # estimate potential outcomes under intervention for CATT x housing vulnerability
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_housing_vulnerability = predict(
      in_housing_spell_on_date_19_model_catt_housing_vulnerability_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x housing vulnerability
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_housing_vulnerability = predict(
      in_housing_spell_on_date_19_model_catt_housing_vulnerability_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 0)
    ),
    # estimate CATT housing vulnerability
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_housing_vulnerability = in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_housing_vulnerability - in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_housing_vulnerability,
    # estimate potential outcomes under intervention for CATT x residential care
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_residential_care = predict(
      in_housing_spell_on_date_19_model_catt_residential_care_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x residential care
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_residential_care = predict(
      in_housing_spell_on_date_19_model_catt_residential_care_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 0)
    ),
    # estimate CATT x residential care
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_residential_care = in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_residential_care - in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_residential_care,
    # estimate potential outcomes under intervention for CATT x placement instability
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_placement_instability = predict(
      in_housing_spell_on_date_19_model_catt_placement_instability_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x placement instability
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_placement_instability = predict(
      in_housing_spell_on_date_19_model_catt_placement_instability_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 0)
    ),
    # estimate CATT x placement instability
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_placement_instability = in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_placement_instability - in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_placement_instability,
    # estimate potential outcomes under intervention for CATT x kinship care
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_kinship_care = predict(
      in_housing_spell_on_date_19_model_catt_kinship_care_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x kinship care
    in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_kinship_care = predict(
      in_housing_spell_on_date_19_model_catt_kinship_care_nn_match,
      type = "response",
      newdata = transform(
        modelling_data_nn_match,
        pyi_flag = 0)
    ),
    # estimate CATT x kinship care
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_kinship_care = in_housing_spell_on_date_19_model_nn_match_potential_outcomes_intervention_kinship_care - in_housing_spell_on_date_19_model_nn_match_potential_outcomes_comparison_kinship_care,
  ) |>
  # subset data
  select(
    pyi_flag,
    distance,
    in_housing_spell_on_date_19_model_nn_match_treatment_effect,
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_sex,
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_aboriginal,
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_housing_vulnerability,
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_residential_care,
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_placement_instability,
    in_housing_spell_on_date_19_model_nn_match_treatment_effect_kinship_care
  ) |>
  pivot_longer(
    cols = c(
      "in_housing_spell_on_date_19_model_nn_match_treatment_effect",
      "in_housing_spell_on_date_19_model_nn_match_treatment_effect_sex",
      "in_housing_spell_on_date_19_model_nn_match_treatment_effect_aboriginal",
      "in_housing_spell_on_date_19_model_nn_match_treatment_effect_housing_vulnerability",
      "in_housing_spell_on_date_19_model_nn_match_treatment_effect_residential_care",
      "in_housing_spell_on_date_19_model_nn_match_treatment_effect_placement_instability",
      "in_housing_spell_on_date_19_model_nn_match_treatment_effect_kinship_care"
    ),
    names_to = "effect_type",
    values_to = "treatment_effect"
  ) |>
  mutate(
    effect_type = case_when(
      effect_type == "in_housing_spell_on_date_19_model_nn_match_treatment_effect" ~ "ATT",
      effect_type == "in_housing_spell_on_date_19_model_nn_match_treatment_effect_sex" ~ "CATT x Sex",
      effect_type == "in_housing_spell_on_date_19_model_nn_match_treatment_effect_aboriginal" ~ "CATT x Aboriginal",
      effect_type == "in_housing_spell_on_date_19_model_nn_match_treatment_effect_housing_vulnerability" ~ "CATT x Housing Vulnerability",
      effect_type == "in_housing_spell_on_date_19_model_nn_match_treatment_effect_residential_care" ~ "CATT x Last placement: residential care",
      effect_type == "in_housing_spell_on_date_19_model_nn_match_treatment_effect_placement_instability" ~ "CATT x History of placement instability",
      effect_type == "in_housing_spell_on_date_19_model_nn_match_treatment_effect_kinship_care" ~ "CATT x Last placement: kinship care"
    ),
    matching_specification = "Nearest neighbour matching"
  )
  
#-------------------------------------------------------------------------------
# 2. Full matching
#-------------------------------------------------------------------------------

# Step A: Fit models
#-----------------------

# ATT model
in_housing_spell_on_date_19_model_att_full_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag  + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_full_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x sex
in_housing_spell_on_date_19_model_catt_sex_full_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * male + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_full_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x aboriginal status
in_housing_spell_on_date_19_model_catt_aboriginal_full_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * aboriginal + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_full_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x housing vulnerability
in_housing_spell_on_date_19_model_catt_housing_vulnerability_full_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * any_housing_spell_before_18 + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_full_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x residential care
in_housing_spell_on_date_19_model_catt_residential_care_full_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * eligible_residential_care_placement + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_full_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x placement instability
in_housing_spell_on_date_19_model_catt_placement_instability_full_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * eligible_placement_instability + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_full_match,
  family = quasibinomial(),
  weights = weights
)

# CATT model x kinship care
in_housing_spell_on_date_19_model_catt_kinship_care_full_match <- glm(
  in_housing_spell_on_date_19 ~ 
    pyi_flag * eligible_kinship_care_placement + 
    eligible_residential_care_placement +
    eligible_time_in_care +                                
    eligible_placement_instability +
    eligible_permanent_placement,
  data = modelling_data_full_match,
  family = quasibinomial(),
  weights = weights
)

# Step B: Estimate Individual TE using G-computation
#---------------------------------------------------

treatment_effect_heterogeneity_in_housing_spell_on_date_19_full_match <- modelling_data_full_match |>
  mutate(
    # estimate potential outcomes under intervention for ATT
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention = predict(
      in_housing_spell_on_date_19_model_att_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for ATT
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison = predict(
      in_housing_spell_on_date_19_model_att_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 0)
    ),
    # estimate ATT
    in_housing_spell_on_date_19_model_full_match_treatment_effect = in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention - in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison,
    # estimate potential outcomes under intervention for CATT x sex
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_sex = predict(
      in_housing_spell_on_date_19_model_catt_sex_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT sex
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_sex = predict(
      in_housing_spell_on_date_19_model_catt_sex_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 0)
    ),
    # estimate CATT sex
    in_housing_spell_on_date_19_model_full_match_treatment_effect_sex = in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_sex - in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_sex,
    # estimate potential outcomes under intervention for CATT x aboriginal status
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_aboriginal = predict(
      in_housing_spell_on_date_19_model_catt_sex_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x aboriginal status
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_aboriginal = predict(
      in_housing_spell_on_date_19_model_catt_sex_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 0)
    ),
    # estimate CATT sex
    in_housing_spell_on_date_19_model_full_match_treatment_effect_aboriginal = in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_aboriginal - in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_aboriginal,
    # estimate potential outcomes under intervention for CATT x housing vulnerability
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_housing_vulnerability = predict(
      in_housing_spell_on_date_19_model_catt_housing_vulnerability_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x housing vulnerability
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_housing_vulnerability = predict(
      in_housing_spell_on_date_19_model_catt_housing_vulnerability_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 0)
    ),
    # estimate CATT housing vulnerability
    in_housing_spell_on_date_19_model_full_match_treatment_effect_housing_vulnerability = in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_housing_vulnerability - in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_housing_vulnerability,
    # estimate potential outcomes under intervention for CATT x residential care
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_residential_care = predict(
      in_housing_spell_on_date_19_model_catt_residential_care_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x residential care
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_residential_care = predict(
      in_housing_spell_on_date_19_model_catt_residential_care_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 0)
    ),
    # estimate CATT x residential care
    in_housing_spell_on_date_19_model_full_match_treatment_effect_residential_care = in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_residential_care - in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_residential_care,
    # estimate potential outcomes under intervention for CATT x placement instability
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_placement_instability = predict(
      in_housing_spell_on_date_19_model_catt_placement_instability_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x placement instability
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_placement_instability = predict(
      in_housing_spell_on_date_19_model_catt_placement_instability_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 0)
    ),
    # estimate CATT x placement instability
    in_housing_spell_on_date_19_model_full_match_treatment_effect_placement_instability = in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_placement_instability - in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_placement_instability,
    # estimate potential outcomes under intervention for CATT x kinship care
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_kinship_care = predict(
      in_housing_spell_on_date_19_model_catt_kinship_care_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 1)
    ),
    # estimate potential outcomes under comparison for CATT x kinship care
    in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_kinship_care = predict(
      in_housing_spell_on_date_19_model_catt_kinship_care_full_match,
      type = "response",
      newdata = transform(
        modelling_data_full_match,
        pyi_flag = 0)
    ),
    # estimate CATT x kinship care
    in_housing_spell_on_date_19_model_full_match_treatment_effect_kinship_care = in_housing_spell_on_date_19_model_full_match_potential_outcomes_intervention_kinship_care - in_housing_spell_on_date_19_model_full_match_potential_outcomes_comparison_kinship_care,
    ) |>
  # subset data
  select(
    pyi_flag,
    distance,
    in_housing_spell_on_date_19_model_full_match_treatment_effect,
    in_housing_spell_on_date_19_model_full_match_treatment_effect_sex,
    in_housing_spell_on_date_19_model_full_match_treatment_effect_aboriginal,
    in_housing_spell_on_date_19_model_full_match_treatment_effect_housing_vulnerability,
    in_housing_spell_on_date_19_model_full_match_treatment_effect_residential_care,
    in_housing_spell_on_date_19_model_full_match_treatment_effect_placement_instability,
    in_housing_spell_on_date_19_model_full_match_treatment_effect_kinship_care
  ) |>
  pivot_longer(
    cols = c(
      "in_housing_spell_on_date_19_model_full_match_treatment_effect",
      "in_housing_spell_on_date_19_model_full_match_treatment_effect_sex",
      "in_housing_spell_on_date_19_model_full_match_treatment_effect_aboriginal",
      "in_housing_spell_on_date_19_model_full_match_treatment_effect_housing_vulnerability",
      "in_housing_spell_on_date_19_model_full_match_treatment_effect_residential_care",
      "in_housing_spell_on_date_19_model_full_match_treatment_effect_placement_instability",
      "in_housing_spell_on_date_19_model_full_match_treatment_effect_kinship_care"
    ),
    names_to = "effect_type",
    values_to = "treatment_effect"
  ) |>
  mutate(
    effect_type = case_when(
      effect_type == "in_housing_spell_on_date_19_model_full_match_treatment_effect" ~ "ATT",
      effect_type == "in_housing_spell_on_date_19_model_full_match_treatment_effect_sex" ~ "CATT x Sex",
      effect_type == "in_housing_spell_on_date_19_model_full_match_treatment_effect_aboriginal" ~ "CATT x Aboriginal",
      effect_type == "in_housing_spell_on_date_19_model_full_match_treatment_effect_housing_vulnerability" ~ "CATT x Housing Vulnerability",
      effect_type == "in_housing_spell_on_date_19_model_full_match_treatment_effect_residential_care" ~ "CATT x Last placement: residential care",
      effect_type == "in_housing_spell_on_date_19_model_full_match_treatment_effect_placement_instability" ~ "CATT x History of placement instability",
      effect_type == "in_housing_spell_on_date_19_model_full_match_treatment_effect_kinship_care" ~ "CATT x Last placement: kinship care"
    ),
    matching_specification = "Full matching"
  )

# Step C: Plot Results
#---------------------------------------------------

#-------------------------------------------------------------------------------
# 2. Combine results for export
#-------------------------------------------------------------------------------

treatment_effect_heterogeneity_in_housing_spell_on_date_19_plot_data <- bind_rows(
  treatment_effect_heterogeneity_in_housing_spell_on_date_19_nn_match,
  treatment_effect_heterogeneity_in_housing_spell_on_date_19_full_match
  )

saveRDS(
  treatment_effect_heterogeneity_in_housing_spell_on_date_19_plot_data,
  "U:/pyi-paper/output/plot_data/hte_diagnostic/hte_in_housing_spell_on_date_19_plot_data.rds"
)
