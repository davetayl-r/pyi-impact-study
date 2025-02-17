# ========================================================================= #
# PYI Impact Study                                                          #
# Matching diagnostic plots                                                 #
# Author: David Taylor                                                      #
# Date: 29/01/2025                                                          #
# ========================================================================= #

# load required packages
library(tidyverse)
library(ggplot2)
library(halfmoon)

# load data
matching_data_file_location <- "P:/pyi/pyi_update/data/processed_data/matching_data.RDS"
matching_data <- readRDS(matching_data_file_location)

full_match_data_location <- "P:/pyi/pyi_update/data/processed_data/pyi_full_match_it_object.RDS"
full_match_data <- readRDS(full_match_data_location)

nn_match_data_location <- "P:/pyi/pyi_update/data/processed_data/pyi_nn_match_it_object.RDS"
nn_match_data <- readRDS(nn_match_data_location)

#---------------------
# Love plot data
#---------------------

# extract info required from object
love_plot_data <- bind_matches(
  matching_data,
  full_match_data,
  nn_match_data
  ) |>
  tidy_smd(
    c(
      male,
      aboriginal,                                               
      eligible_residential_care_placement,                     
      eligible_time_in_care,                      
      eligible_placement_instability,
      eligible_permanent_placement,
      eligible_kinship_care_placement,   
      permanent_responsibility_minister_eligiblity_period,
      any_housing_spell_before_18
    ),
    .group = pyi_flag,
    .wts = c(
      full_match_data,
      nn_match_data
    )
  )

# export love plot data
saveRDS(
  love_plot_data,
  "./output/plot_data/love_plot_data.RDS"
)

#---------------------
# Balance plot data
#---------------------

# nn match
balance_plot_data_nn_match <- bind_matches(
  matching_data,
  nn_match_data
  ) 

balance_plot_data_nn_match$ps_score <- nn_match_data$distance
balance_plot_data_nn_match$weights <- nn_match_data$weights
balance_plot_data_nn_match$matching_specification <- "Nearest Neighbour"

balance_plot_data_nn_match_subset <- balance_plot_data_nn_match |>
  select(
    pyi_flag,
    ps_score,
    weights,
    matching_specification
  )

# full match
balance_plot_data_full_match <- bind_matches(
  matching_data,
  full_match_data
) 

balance_plot_data_full_match$ps_score <- full_match_data$distance
balance_plot_data_full_match$weights <- full_match_data$weights
balance_plot_data_full_match$matching_specification <- "Full"

balance_plot_data_full_match_subset <- balance_plot_data_full_match |>
  select(
    pyi_flag,
    ps_score,
    weights,
    matching_specification
  )

# combine full and nn match
balance_plot_data <- bind_rows(
  balance_plot_data_nn_match_subset,
  balance_plot_data_full_match_subset
)

# export balance plot data
saveRDS(
  balance_plot_data,
  "./output/plot_data/balance_plot_data.RDS"
)

