# ========================================================================= #
# PYI Impact Study                                                          #
# Create matching data set                                                  #
# Author: David Taylor                                                      #
# Date: 15/01/2024                                                          #
# ========================================================================= #

# This code merges summarised oohc and cims files to create and analysis dataset 

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages 
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)

# set up workspace preferences using custom helper function
workspace_setup()

# read oohc data
summarised_oohc_file_location <- "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_fifteen.RDS"
summarised_oohc_file <- readRDS(summarised_oohc_file_location)

# read cims data
summarised_cims_file_location <- "P:/pyi/pyi_update/data/processed_data/cims_data_summarise_step_six.RDS"
summarised_cims_file <- readRDS(summarised_cims_file_location)

# merge files
merged_data <- left_join(
  summarised_oohc_file,
  summarised_cims_file,
  by = "childstory_id"
  )

# clear up data
clean_merged_data <- merged_data |>
  mutate(
    # replace na values
    across(
      .cols = c(
        "any_housing_spell_before_18",
        "in_housing_spell_on_date_18",
        "in_housing_spell_on_date_19",
        "new_housing_spell_between_18_19",
        "any_housing_spell_between_18_19",                      
        "time_in_housing_support_18_19",
        "new_unsheltered_homelessness_between_18_19",
        "new_ongoing_unsheltered_homelessness_between_18_19",
        "new_short_term_accommodation_required_between_18_19",    
        "new_ongoing_short_term_accommodation_required_between_18_19",
        "new_requires_housing_assistance_between_18_19",
        "new_ongoing_requires_housing_assistance_between_18_19",
        "count_total_spells_18_19"),
      .fns = ~replace_na(.,0)
    )
  ) |>
  # drop redundant vars
  select(
    -aihw_slk.x,
    -aihw_slk.y,
    -slk_trimmed,
    -district_catchment,
    -pyi_eligible_location,
    -in_pyi_location_did_not_receive,
    -matching_pool_flag,
    -date_19,
    -date_18,
    -date_16,
    -residential_care_placement,
    -eligible_recent_placement_instability,
    -eligible_lifetime_placement_instability,
    -self_placed_missing_absent_from_placement,
    -independent_living_placement,
    -eligibility_window_start,
    -eligibility_window_end,
    -days_eligible_pyi,
    -meets_pyi_eligibility_criteria,
    -meets_pyi_eligibility_criteria_during_evaluation,
    -placement_ended_disruptive_behaviour,
    -placement_breakdown,
    -placement_breakdown_allegation_against_carer
  ) |>
  ungroup()

# export analysis data
saveRDS(
  clean_merged_data, 
  "P:/pyi/pyi_update/data/processed_data/matching_data.RDS")