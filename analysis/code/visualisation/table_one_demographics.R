# ========================================================================= #
# PYI Impact Study                                                          #
# Table One                                                                 #
# Author: David Taylor                                                      #
# Date: 29/01/2025                                                          #
# ========================================================================= #

# load required packages
library(tidyverse)

# read data
matching_data_file_location <- "P:/pyi/pyi_update/data/processed_data/matching_data.RDS"
matching_data <- readRDS(matching_data_file_location)

# create data for table one
table_one_data_binary <- matching_data |>
  # select only pyi sample
  filter(
    pyi_flag == 1
  ) |>
  # select vars for table
  select(
    male,                                                        
    female,                                                     
    aboriginal,
    eligible_residential_care_placement,                         
    eligible_foster_care_placement,                          
    eligible_kinship_care_placement,                             
    eligible_permanent_placement,                              
    permanent_responsibility_minister_eligiblity_period,         
    independent_living_placement_after_16,        
    self_placed_missing_absent_from_placement_after_16,
    eligible_time_in_care,                          
    eligible_placement_instability,                             
    last_placement_self_placed_missing_absent,                   
    last_placement_independent_living,                    
    placement_ended_disruptive_behaviour_after_16,               
    placement_breakdown_after_16,             
    any_housing_spell_before_18
  ) |>
  mutate(
    non_aboriginal = case_when(
      aboriginal == 0 ~ 1,
      TRUE ~ 0
    )
  ) |>
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) |>
  group_by(
    variable
  ) |>
  summarise(
    n = sum(
      value == 1,
      na.rm = TRUE
      ),
    total = n(),
    percent = round(n/total * 100, 2)
  ) |>
  mutate(
    n_percent = paste(n, " (", percent, ")", sep = "")
  ) |>
  select(
    -n,
    -percent
  )
  
# create data for table one
time_in_care_mean_sd <- matching_data |>
  # select only pyi sample
  filter(
    pyi_flag == 1
  ) |>
  select(
    time_in_care_start_eligible_window                      
  ) |>
  summarise(
    time_in_care_mean = round(
      mean(time_in_care_start_eligible_window), 
      2),
    time_in_care_sd = round(
      sd(time_in_care_start_eligible_window), 
      2)
  )
  
table_one_data_continuous <- tibble(
    variable = "time_in_care_start_eligible_window (mean (SD))",
    mean_sd = paste(time_in_care_mean_sd$time_in_care_mean, " (", time_in_care_mean_sd$time_in_care_sd, ")", sep = "")
  )

# join binary and continuous data together
table_one_data <- full_join(
  table_one_data_binary,
  table_one_data_continuous
  )
  
# export data
write_csv(
  table_one_data,
  "U:/pyi-paper/output/table_data/table_one_data.csv"
)