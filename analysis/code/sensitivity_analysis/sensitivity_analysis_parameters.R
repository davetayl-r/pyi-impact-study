# parameters for sensitivity analysis

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages using custom helper function
library(tidyverse)
library(MatchIt)
library(cobalt)
library(tableone)
library(tinytable)

# set up workspace preferences using custom helper function
workspace_setup()

# read data
matching_data_file_location <- "P:/pyi/pyi_update/data/processed_data/matching_data.RDS"
matching_data <- readRDS(matching_data_file_location)

# prepare data for table one
table_one_unmatched_data <- matching_data |>
  select(
    pyi_flag,
    male,
    aboriginal,                                                
    eligible_residential_care_placement,                       
    eligible_time_in_care,                                   
    eligible_placement_instability,
    eligible_permanent_placement,
    eligible_kinship_care_placement,       
    permanent_responsibility_minister_eligiblity_period,
    time_in_care_start_eligible_window,
    last_placement_self_placed_missing_absent,                 
    last_placement_independent_living,              
    self_placed_missing_absent_from_placement_after_16,         
    independent_living_placement_after_16,
    any_housing_spell_before_18
  ) |>
  mutate(
    across(
      c(
        pyi_flag,
        male,
        aboriginal,                                                
        eligible_residential_care_placement,                       
        eligible_time_in_care,                                   
        eligible_placement_instability,
        eligible_permanent_placement,
        eligible_kinship_care_placement,       
        permanent_responsibility_minister_eligiblity_period,
        last_placement_self_placed_missing_absent,                 
        last_placement_independent_living,              
        self_placed_missing_absent_from_placement_after_16,         
        independent_living_placement_after_16,
        any_housing_spell_before_18
      ),
      as.factor)
  )

# get list of vars to show in table one
variable_list <- c(
  "male",
  "aboriginal",                                                
  "eligible_residential_care_placement",                       
  "eligible_time_in_care",                                   
  "eligible_placement_instability",
  "eligible_permanent_placement",
  "eligible_kinship_care_placement",       
  "permanent_responsibility_minister_eligiblity_period",
  "any_housing_spell_before_18"
)

# create table one for unmatched data
table_one_unmatched <- CreateTableOne(
  vars = variable_list,
  strata = "pyi_flag",
  test = FALSE,
  smd = TRUE,
  data = table_one_unmatched_data
)

# print and convert to data.frame
table_one_unmatched_data_output <- print(
  table_one_unmatched,
  smd = TRUE) |>
  as.data.frame() |>
  rownames_to_column(
    var = "Variable"
  ) |>
  rename(
    `Comparison Pool` = 2,
    `Received Intervention` = 3,
    `SMD\n(unmatched)` = 4
  )

# get minimum and maximum prevalence from adjustment set
tip_analysis_inputs <- table_one_unmatched_data_output |>
  filter(
    !Variable == "n"
  ) |>
  select(
    Variable,
    `Received Intervention`,
    `Comparison Pool`
  ) |>
  rename(
    received_intervention = "Received Intervention",
    comparison_pool = "Comparison Pool"
  ) |>
  mutate(
    received_intervention_percent = str_extract(
      received_intervention, "\\((.*?)\\)"),
    received_intervention_percent = str_remove_all(
      received_intervention_percent, "[()]"),
    received_intervention_percent = as.numeric(received_intervention_percent),
    comparison_pool_percent = str_extract(
      comparison_pool, "\\((.*?)\\)"),
    comparison_pool_percent = str_remove_all(
      comparison_pool_percent, "[()]"),
    comparison_pool_percent = as.numeric(comparison_pool_percent),
    prevalence_treated = received_intervention_percent/100,
    prevalence_untreated = comparison_pool_percent/100,
    difference_prevalence = abs(prevalence_untreated - prevalence_treated)
  ) |>
  select(
    prevalence_treated,
    prevalence_untreated,
    difference_prevalence
  ) |>
  arrange(
    -prevalence_treated
  ) |>
  slice(
    c(3, 4, 9)
  ) |>
  mutate(
    # formulate into scenarios
    scenario = case_when(
      row_number() == 1 ~ "high_prevalence",
      row_number() == 2 ~ "medium_prevalence",
      row_number() == 3 ~ "low_prevalence"
    ),
    # round prevalence_treated values
    prevalence_treated = case_when(
      scenario == "high_prevalence" ~ 0.95,
      scenario == "medium_prevalence" ~ 0.75,
      scenario == "low_prevalence" ~ 0.25,
    ),
    # use max difference_prevalence for all scenarios
    difference_prevalence = 0.2,
    # calculate difference_prevalence
    prevalence_untreated = prevalence_treated - difference_prevalence
  ) |>
  # reorder data.frame
  select(
    scenario,
    prevalence_treated,
    prevalence_untreated,
    difference_prevalence
  )

# export for sensitivity analysis
saveRDS(
  tip_analysis_inputs,
  "U:/pyi-paper/output/tip_analysis_inputs.RDS"
)

