# ========================================================================= #
# PYI Impact Study                                                          #
# Run nearest neighbour matching model                                      #
# Author: David Taylor                                                      #
# Date: 29/01/2025                                                          #
# ========================================================================= #

# This code:
#   1. Fits a nearest neighbour matching model in MatchIt
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
# 1. Run NN Matching Model using MatchIt 
#-------------------------------------------------------------------------------

# run matching model informed by dag
pyi_dag_match_nn <- matchit(
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
  method = "nearest",
  distance = "gam",
  link = "probit",
  data = matching_data,
  estimand = "ATT"
)

#-------------------------------------------------------------------------------
# 2. Inspect results 
#-------------------------------------------------------------------------------

# inspect results
pyi_dag_match_nn_summary <- summary(
  pyi_dag_match_nn,
  un = FALSE
  )
pyi_dag_match_nn_summary

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

# inspect love plot from cobalt
cobalt::love.plot(
  pyi_dag_match_nn,
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

# inspect balance plot from cobalt
cobalt::bal.plot(
  pyi_dag_match_nn,
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

saveRDS(
  pyi_dag_match_nn,
  "P:/pyi/pyi_update/data/processed_data/pyi_nn_match_it_object.RDS"
)

#-------------------------------------------------------------------------------
# 4. Extract matched data and export results 
#-------------------------------------------------------------------------------

# extract matched data
nn_matched_data <- match.data(
  pyi_dag_match_nn
)

# export matched data for modelling
saveRDS(
  nn_matched_data, 
  "P:/pyi/pyi_update/data/processed_data/pyi_modelling_data_nn_match.RDS"
  )
