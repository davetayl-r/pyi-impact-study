# ========================================================================= #
# PYI Impact Study                                                          #
# Run full matching model                                                   #
# Author: David Taylor                                                      #
# Date: 29/01/2025                                                          #
# ========================================================================= #

# This code:
#   1. Fits a full matching model in MatchIt
#   2. Exports the results 

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages using custom helper function
library(tidyverse)
library(MatchIt)
library(cobalt)

# set up workspace preferences using custom helper function
workspace_setup()

# read data
matching_data_file_location <- "P:/pyi/pyi_update/data/processed_data/matching_data.RDS"
matching_data <- readRDS(matching_data_file_location)

#-------------------------------------------------------------------------------
# 1. Run Full Matching Model using MatchIt 
#-------------------------------------------------------------------------------

# run full matching model informed by dag
pyi_dag_full_match <- matchit(
  pyi_flag ~ 
    male +
    aboriginal +                                                
    eligible_residential_care_placement +                       
    eligible_time_in_care +                                   
    eligible_placement_instability +
    eligible_permanent_placement +
    eligible_kinship_care_placement +       
    permanent_responsibility_minister_eligiblity_period +
    any_housing_spell_before_18,
  replace = FALSE,
  method = "full",
  discard = "control",
  distance = "gam",
  link = "probit",
  data = matching_data,
  estimand = "ATT"
)

#-------------------------------------------------------------------------------
# 2. Inspect results 
#-------------------------------------------------------------------------------

# inspect results
pyi_dag_full_match_summary <- summary(
  pyi_dag_full_match,
  un = FALSE
  )
pyi_dag_full_match_summary

# rename variable names
variable_names <- data.frame(
  old = c(
    "distance",
    "male",
    "aboriginal",                        
    "eligible_residential_care_placement",                       
    "eligible_time_in_care",                           
    "eligible_placement_instability",
    "eligible_permanent_placement",
    "eligible_kinship_care_placement",
    "permanent_responsibility_minister_eligiblity_period",
    "last_placement_self_placed_missing_absent",                
    "any_housing_spell_before_18"
  ),
  new = c(
    "Propensity Score",
    "Male",
    "Aboriginal",                                                
    "In residential care placement during eligiblity period",                       
    "Twelve months or more in OOHC by end of eligibility period",                                   
    "Current or previous placement instability",
    "In permanent OOHC placement during eligiblity period",
    "In kinship care placement during eligiblity period",       
    "Permanent Responsibility of the Minister during eligiblity period",
    "Last placement before 18 was self-placed, missing or absent",                 
    "Use of specialist homelessness services between age 16 and 18"
  )
)

# inspect love plot
cobalt::love.plot(
  pyi_dag_full_match,
  stats = c("m"),
  thresholds = c(m = 0.1),
  binary = "std",
  var.names = variable_names,
  var.order = "unadjusted",
  colours = c(
    "#EE4266",
    "#0EAD69"
  ),
  shapes = c(
    "square",
    "triangle"
  ),
  abs = TRUE,
  position = "bottom",
  xlab = "Absolute Standardised Mean Difference"
  )

# inspect balance plot
cobalt::bal.plot(
  pyi_dag_full_match,
  which = "both",
  var.name = "distance",
  type = "histogram",
  mirror = TRUE,
  colours = c(
    "#EE4266",
    "#0EAD69"
  )
  )

#-------------------------------------------------------------------------------
# 3. Export MatchIt object for plotting  
#-------------------------------------------------------------------------------

# save matchit object for plotting
saveRDS(
  pyi_dag_full_match,
  "P:/pyi/pyi_update/data/processed_data/pyi_full_match_it_object.RDS"
)

#-------------------------------------------------------------------------------
# 4. Extract matched data and export results 
#-------------------------------------------------------------------------------

# extract matched data
full_matched_data <- match.data(
  pyi_dag_full_match
)

# export matched data for modelling
saveRDS(
  full_matched_data, 
  "P:/pyi/pyi_update/data/processed_data/pyi_modelling_data_full_match.RDS"
  )
