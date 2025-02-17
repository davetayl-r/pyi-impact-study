# ========================================================================= #
# PYI Impact Study                                                          #
# Visualisation data preparation                                             #
# Author: David Taylor                                                      #
# Date: 23/01/2025                                                          #
# ========================================================================= #

# This code prepares output from this analysis for export, this includes (for both nn and full match specifications): 
#   1. Treatment effect results 
#   2. Subgroup analysis x sex results
#   3. Subgroup analysis x aboriginal results
#   4. Subgroup analysis x vulnerability results
#   5. Sensitivity analysis results 
#   6. SMD Summmary plot

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages
library(tidyverse)

# set up workspace preferences using custom helper function
workspace_setup()

#-------------------------------------------------------------------------------
# 2. Table 2: Treatment effect results
#-------------------------------------------------------------------------------

# read in all treatment effect results
sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "te_results_nn_match"
)

sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "te_results_full_match"
)

# combine results from all models
combined_att_results <- bind_cols(
    in_housing_spell_on_date_18_te_results_full_match,
    in_housing_spell_on_date_18_te_results_nn_match,
    in_housing_spell_on_date_19_te_results_full_match,
    in_housing_spell_on_date_19_te_results_nn_match,
    new_housing_spell_between_18_19_te_results_full_match,
    new_housing_spell_between_18_19_te_results_nn_match,
    any_housing_spell_between_18_19_te_results_full_match,
    any_housing_spell_between_18_19_te_results_nn_match,
    new_unsheltered_homelessness_between_18_19_te_results_full_match,
    new_unsheltered_homelessness_between_18_19_te_results_nn_match,
    new_ongoing_unsheltered_homelessness_between_18_19_te_results_full_match,
    new_ongoing_unsheltered_homelessness_between_18_19_te_results_nn_match,
    new_short_term_accommodation_required_between_18_19_te_results_full_match,
    new_short_term_accommodation_required_between_18_19_te_results_nn_match,
    new_ongoing_short_term_accommodation_required_between_18_19_te_results_full_match,
    new_ongoing_short_term_accommodation_required_between_18_19_te_results_nn_match,
    count_total_spells_18_19_te_results_full_match,
    count_total_spells_18_19_te_results_nn_match,
    time_in_housing_support_18_19_te_results_full_match,
    time_in_housing_support_18_19_te_results_nn_match
  ) |>
  rename(
    metrics = 1,
    in_housing_spell_on_date_18_full_match = in_housing_spell_on_date_18...2,                                 
    in_housing_spell_on_date_18_nn_match = in_housing_spell_on_date_18...4,                                 
    in_housing_spell_on_date_19_full_match = in_housing_spell_on_date_19...6,                                 
    in_housing_spell_on_date_19_nn_match = in_housing_spell_on_date_19...8,                            
    new_housing_spell_between_18_19_full_match = new_housing_spell_between_18_19...10,                            
    new_housing_spell_between_18_19_nn_match = new_housing_spell_between_18_19...12,                            
    any_housing_spell_between_18_19_full_match = any_housing_spell_between_18_19...14,                            
    any_housing_spell_between_18_19_nn_match = any_housing_spell_between_18_19...16,                 
    new_unsheltered_homelessness_between_18_19_full_match = new_unsheltered_homelessness_between_18_19...18,                 
    new_unsheltered_homelessness_between_18_19_nn_match = new_unsheltered_homelessness_between_18_19...20,         
    new_ongoing_unsheltered_homelessness_between_18_19_full_match = new_ongoing_unsheltered_homelessness_between_18_19...22,         
    new_ongoing_unsheltered_homelessness_between_18_19_nn_match = new_ongoing_unsheltered_homelessness_between_18_19...24,          
    new_short_term_accommodation_required_between_18_19_full_match = new_short_term_accommodation_required_between_18_19...26,        
    new_short_term_accommodation_required_between_18_19_nn_match = new_short_term_accommodation_required_between_18_19...28,
    new_ongoing_short_term_accommodation_required_between_18_19_full_match = new_ongoing_short_term_accommodation_required_between_18_19...30,
    new_ongoing_short_term_accommodation_required_between_18_19_nn_match = new_ongoing_short_term_accommodation_required_between_18_19...32,
    count_total_spells_18_19_full_match = count_total_spells_18_19...34,                                   
    count_total_spells_18_19_nn_match = count_total_spells_18_19...36,                              
    time_in_housing_support_18_19_full_match = time_in_housing_support_18_19...38,                              
    time_in_housing_support_18_19_nn_match = time_in_housing_support_18_19...40
  ) |>
  select(
    -starts_with("reported_metrics...")
  )

# export results
write_csv(
  combined_att_results,
  "U:/pyi-paper/output/treatment_effect_results/combined_att_results.csv"
)

#-------------------------------------------------------------------------------
# 3. Table 3: Conditional treatment effect results: x Sex
#-------------------------------------------------------------------------------

# read in all subgroup analysis x sex results
sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "subgroup_sex_results_full_match"
)

sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "subgroup_sex_results_nn_match"
)

# combine results from all models
combined_subgroup_sex_results <- bind_cols(
    in_housing_spell_on_date_18_subgroup_sex_results_full_match,
    in_housing_spell_on_date_18_subgroup_sex_results_nn_match,
    in_housing_spell_on_date_19_subgroup_sex_results_full_match,
    in_housing_spell_on_date_19_subgroup_sex_results_nn_match,
    new_housing_spell_between_18_19_subgroup_sex_results_full_match,
    new_housing_spell_between_18_19_subgroup_sex_results_nn_match,
    any_housing_spell_between_18_19_subgroup_sex_results_full_match,
    any_housing_spell_between_18_19_subgroup_sex_results_nn_match,
    new_unsheltered_homelessness_between_18_19_subgroup_sex_results_full_match,
    new_unsheltered_homelessness_between_18_19_subgroup_sex_results_nn_match,
    new_ongoing_unsheltered_homelessness_between_18_19_subgroup_sex_results_full_match,
    new_ongoing_unsheltered_homelessness_between_18_19_subgroup_sex_results_nn_match,
    new_short_term_accommodation_required_between_18_19_subgroup_sex_results_full_match,
    new_short_term_accommodation_required_between_18_19_subgroup_sex_results_nn_match,
    new_ongoing_short_term_accommodation_required_between_18_19_subgroup_sex_results_full_match,
    new_ongoing_short_term_accommodation_required_between_18_19_subgroup_sex_results_nn_match,
    count_total_spells_18_19_subgroup_sex_results_full_match,
    count_total_spells_18_19_subgroup_sex_results_nn_match,
    time_in_housing_support_18_19_subgroup_sex_results_full_match,
    time_in_housing_support_18_19_subgroup_sex_results_nn_match
  ) |>
  rename(
    metrics = 1,
    in_housing_spell_on_date_18_male_full_match = in_housing_spell_on_date_18_male...2,                             
    in_housing_spell_on_date_18_female_full_match = in_housing_spell_on_date_18_female...3,                                
    in_housing_spell_on_date_18_male_nn_match = in_housing_spell_on_date_18_male...5,                                  
    in_housing_spell_on_date_18_female_nn_match = in_housing_spell_on_date_18_female...6,                                 
    in_housing_spell_on_date_19_male_full_match = in_housing_spell_on_date_19_male...8,                                 
    in_housing_spell_on_date_19_female_full_match = in_housing_spell_on_date_19_female...9,                                  
    in_housing_spell_on_date_19_male_nn_match = in_housing_spell_on_date_19_male...11,                                 
    in_housing_spell_on_date_19_female_nn_match = in_housing_spell_on_date_19_female...12,                                 
    new_housing_spell_between_18_19_male_full_match = new_housing_spell_between_18_19_male...14,                      
    new_housing_spell_between_18_19_female_full_match = new_housing_spell_between_18_19_female...15,                             
    new_housing_spell_between_18_19_male_nn_match = new_housing_spell_between_18_19_male...17,                       
    new_housing_spell_between_18_19_female_nn_match = new_housing_spell_between_18_19_female...18,                             
    any_housing_spell_between_18_19_male_full_match = any_housing_spell_between_18_19_male...20,                       
    any_housing_spell_between_18_19_female_full_match = any_housing_spell_between_18_19_female...21,                             
    any_housing_spell_between_18_19_male_nn_match = any_housing_spell_between_18_19_male...23,                       
    any_housing_spell_between_18_19_female_nn_match = any_housing_spell_between_18_19_female...24,                         
    new_unsheltered_homelessness_between_18_19_male_full_match = new_unsheltered_homelessness_between_18_19_male...26,            
    new_unsheltered_homelessness_between_18_19_female_full_match = new_unsheltered_homelessness_between_18_19_female...27,                  
    new_unsheltered_homelessness_between_18_19_male_nn_match = new_unsheltered_homelessness_between_18_19_male...29,            
    new_unsheltered_homelessness_between_18_19_female_nn_match = new_unsheltered_homelessness_between_18_19_female...30,                 
    new_ongoing_unsheltered_homelessness_between_18_19_male_full_match = new_ongoing_unsheltered_homelessness_between_18_19_male...32,    
    new_ongoing_unsheltered_homelessness_between_18_19_female_full_match = new_ongoing_unsheltered_homelessness_between_18_19_female...33,          
    new_ongoing_unsheltered_homelessness_between_18_19_male_nn_match = new_ongoing_unsheltered_homelessness_between_18_19_male...35,    
    new_ongoing_unsheltered_homelessness_between_18_19_female_nn_match = new_ongoing_unsheltered_homelessness_between_18_19_female...36,          
    new_short_term_accommodation_required_between_18_19_male_full_match = new_short_term_accommodation_required_between_18_19_male...38,   
    new_short_term_accommodation_required_between_18_19_female_full_match = new_short_term_accommodation_required_between_18_19_female...39,         
    new_short_term_accommodation_required_between_18_19_male_nn_match = new_short_term_accommodation_required_between_18_19_male...41,   
    new_short_term_accommodation_required_between_18_19_female_nn_match = new_short_term_accommodation_required_between_18_19_female...42,    
    new_ongoing_short_term_accommodation_required_between_18_19_male_full_match = new_ongoing_short_term_accommodation_required_between_18_19_male...44,
    new_ongoing_short_term_accommodation_required_between_18_19_female_full_match = new_ongoing_short_term_accommodation_required_between_18_19_female...45, 
    new_ongoing_short_term_accommodation_required_between_18_19_male_nn_match = new_ongoing_short_term_accommodation_required_between_18_19_male...47,
    new_ongoing_short_term_accommodation_required_between_18_19_female_nn_match = new_ongoing_short_term_accommodation_required_between_18_19_female...48, 
    count_total_spells_18_19_male_full_match = count_total_spells_18_19_male...50,                              
    count_total_spells_18_19_female_full_match = count_total_spells_18_19_female...51,                                    
    count_total_spells_18_19_male_nn_match = count_total_spells_18_19_male...53,                              
    count_total_spells_18_19_female_nn_match = count_total_spells_18_19_female...54,                                    
    time_in_housing_support_18_19_male_full_match = time_in_housing_support_18_19_male...56,                           
    time_in_housing_support_18_19_female_full_match = time_in_housing_support_18_19_female...57,                               
    time_in_housing_support_18_19_male_nn_match = time_in_housing_support_18_19_male...59,                          
    time_in_housing_support_18_19_female_nn_match = time_in_housing_support_18_19_female...60 
  ) |>
  select(
    -starts_with("reported_metrics...")
  )

# export results
write_csv(
  combined_subgroup_sex_results,
  "U:/pyi-paper/output/treatment_effect_results/combined_subgroup_sex_results.csv"
)

#-------------------------------------------------------------------------------
# 4. Table 4: Conditional treatment effect results: x Aboriginal status
#-------------------------------------------------------------------------------

# read in all subgroup analysis x sex results
sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "subgroup_aboriginal_results_full_match"
)

sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "subgroup_aboriginal_results_nn_match"
)

# combine results from all models
combined_subgroup_aboriginal_results <- bind_cols(
  in_housing_spell_on_date_18_subgroup_aboriginal_results_full_match,
  in_housing_spell_on_date_18_subgroup_aboriginal_results_nn_match,
  in_housing_spell_on_date_19_subgroup_aboriginal_results_full_match,
  in_housing_spell_on_date_19_subgroup_aboriginal_results_nn_match,
  new_housing_spell_between_18_19_subgroup_aboriginal_results_full_match,
  new_housing_spell_between_18_19_subgroup_aboriginal_results_nn_match,
  any_housing_spell_between_18_19_subgroup_aboriginal_results_full_match,
  any_housing_spell_between_18_19_subgroup_aboriginal_results_nn_match,
  new_unsheltered_homelessness_between_18_19_subgroup_aboriginal_results_full_match,
  new_unsheltered_homelessness_between_18_19_subgroup_aboriginal_results_nn_match,
  new_ongoing_unsheltered_homelessness_between_18_19_subgroup_aboriginal_results_full_match,
  new_ongoing_unsheltered_homelessness_between_18_19_subgroup_aboriginal_results_nn_match,
  new_short_term_accommodation_required_between_18_19_subgroup_aboriginal_results_full_match,
  new_short_term_accommodation_required_between_18_19_subgroup_aboriginal_results_nn_match,
  new_ongoing_short_term_accommodation_required_between_18_19_subgroup_aboriginal_results_full_match,
  new_ongoing_short_term_accommodation_required_between_18_19_subgroup_aboriginal_results_nn_match,
  count_total_spells_18_19_subgroup_aboriginal_results_full_match,
  count_total_spells_18_19_subgroup_aboriginal_results_nn_match,
  time_in_housing_support_18_19_subgroup_aboriginal_results_full_match,
  time_in_housing_support_18_19_subgroup_aboriginal_results_nn_match
) |>
  rename(
    metrics = 1,
    in_housing_spell_on_date_18_aboriginal_full_match = in_housing_spell_on_date_18_aboriginal...2,                             
    in_housing_spell_on_date_18_non_aboriginal_full_match = in_housing_spell_on_date_18_non_aboriginal...3,                                
    in_housing_spell_on_date_18_aboriginal_nn_match = in_housing_spell_on_date_18_aboriginal...5,                                  
    in_housing_spell_on_date_18_non_aboriginal_nn_match = in_housing_spell_on_date_18_non_aboriginal...6,                                 
    in_housing_spell_on_date_19_aboriginal_full_match = in_housing_spell_on_date_19_aboriginal...8,                                 
    in_housing_spell_on_date_19_non_aboriginal_full_match = in_housing_spell_on_date_19_non_aboriginal...9,                                  
    in_housing_spell_on_date_19_aboriginal_nn_match = in_housing_spell_on_date_19_aboriginal...11,                                 
    in_housing_spell_on_date_19_non_aboriginal_nn_match = in_housing_spell_on_date_19_non_aboriginal...12,                                 
    new_housing_spell_between_18_19_aboriginal_full_match = new_housing_spell_between_18_19_aboriginal...14,                      
    new_housing_spell_between_18_19_non_aboriginal_full_match = new_housing_spell_between_18_19_non_aboriginal...15,                             
    new_housing_spell_between_18_19_aboriginal_nn_match = new_housing_spell_between_18_19_aboriginal...17,                       
    new_housing_spell_between_18_19_non_aboriginal_nn_match = new_housing_spell_between_18_19_non_aboriginal...18,                             
    any_housing_spell_between_18_19_aboriginal_full_match = any_housing_spell_between_18_19_aboriginal...20,                       
    any_housing_spell_between_18_19_non_aboriginal_full_match = any_housing_spell_between_18_19_non_aboriginal...21,                             
    any_housing_spell_between_18_19_aboriginal_nn_match = any_housing_spell_between_18_19_aboriginal...23,                       
    any_housing_spell_between_18_19_non_aboriginal_nn_match = any_housing_spell_between_18_19_non_aboriginal...24,                         
    new_unsheltered_homelessness_between_18_19_aboriginal_full_match = new_unsheltered_homelessness_between_18_19_aboriginal...26,            
    new_unsheltered_homelessness_between_18_19_non_aboriginal_full_match = new_unsheltered_homelessness_between_18_19_non_aboriginal...27,                  
    new_unsheltered_homelessness_between_18_19_aboriginal_nn_match = new_unsheltered_homelessness_between_18_19_aboriginal...29,            
    new_unsheltered_homelessness_between_18_19_non_aboriginal_nn_match = new_unsheltered_homelessness_between_18_19_non_aboriginal...30,                 
    new_ongoing_unsheltered_homelessness_between_18_19_aboriginal_full_match = new_ongoing_unsheltered_homelessness_between_18_19_aboriginal...32,    
    new_ongoing_unsheltered_homelessness_between_18_19_non_aboriginal_full_match = new_ongoing_unsheltered_homelessness_between_18_19_non_aboriginal...33,          
    new_ongoing_unsheltered_homelessness_between_18_19_aboriginal_nn_match = new_ongoing_unsheltered_homelessness_between_18_19_aboriginal...35,    
    new_ongoing_unsheltered_homelessness_between_18_19_non_aboriginal_nn_match = new_ongoing_unsheltered_homelessness_between_18_19_non_aboriginal...36,          
    new_short_term_accommodation_required_between_18_19_aboriginal_full_match = new_short_term_accommodation_required_between_18_19_aboriginal...38,   
    new_short_term_accommodation_required_between_18_19_non_aboriginal_full_match = new_short_term_accommodation_required_between_18_19_non_aboriginal...39,         
    new_short_term_accommodation_required_between_18_19_aboriginal_nn_match = new_short_term_accommodation_required_between_18_19_aboriginal...41,   
    new_short_term_accommodation_required_between_18_19_non_aboriginal_nn_match = new_short_term_accommodation_required_between_18_19_non_aboriginal...42,        
    new_ongoing_short_term_accommodation_required_between_18_19_aboriginal_full_match = new_ongoing_short_term_accommodation_required_between_18_19_aboriginal...44,
    new_ongoing_short_term_accommodation_required_between_18_19_non_aboriginal_full_match = new_ongoing_short_term_accommodation_required_between_18_19_non_aboriginal...45, 
    new_ongoing_short_term_accommodation_required_between_18_19_aboriginal_nn_match = new_ongoing_short_term_accommodation_required_between_18_19_aboriginal...47,
    new_ongoing_short_term_accommodation_required_between_18_19_non_aboriginal_nn_match = new_ongoing_short_term_accommodation_required_between_18_19_non_aboriginal...48, 
    count_total_spells_18_19_aboriginal_full_match = count_total_spells_18_19_aboriginal...50,                              
    count_total_spells_18_19_non_aboriginal_full_match = count_total_spells_18_19_non_aboriginal...51,                                    
    count_total_spells_18_19_aboriginal_nn_match = count_total_spells_18_19_aboriginal...53,                              
    count_total_spells_18_19_non_aboriginal_nn_match = count_total_spells_18_19_non_aboriginal...54,                                    
    time_in_housing_support_18_19_aboriginal_full_match = time_in_housing_support_18_19_aboriginal...56,                           
    time_in_housing_support_18_19_non_aboriginal_full_match = time_in_housing_support_18_19_non_aboriginal...57,                               
    time_in_housing_support_18_19_aboriginal_nn_match = time_in_housing_support_18_19_aboriginal...59,                          
    time_in_housing_support_18_19_non_aboriginal_nn_match = time_in_housing_support_18_19_non_aboriginal...60 
  ) |>
  select(
    -starts_with("reported_metrics...")
  )

# export results
write_csv(
  combined_subgroup_aboriginal_results,
  "U:/pyi-paper/output/treatment_effect_results/combined_subgroup_aboriginal_results.csv"
)

#-------------------------------------------------------------------------------
# 5. Table 5: Conditional treatment effect results: x Housing vulnerability
#-------------------------------------------------------------------------------

# read in all subgroup analysis x housing vulnerability results
sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "subgroup_housing_vulnerability_results_full_match"
)

sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/treatment_effect_results/",
  common_suffix = "subgroup_housing_vulnerability_results_nn_match"
)

# combine results from all models
combined_subgroup_housing_vulnerability_results <- bind_cols(
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_results_full_match,
  in_housing_spell_on_date_18_subgroup_housing_vulnerability_results_nn_match,
  in_housing_spell_on_date_19_subgroup_housing_vulnerability_results_full_match,
  in_housing_spell_on_date_19_subgroup_housing_vulnerability_results_nn_match,
  new_housing_spell_between_18_19_subgroup_housing_vulnerability_results_full_match,
  new_housing_spell_between_18_19_subgroup_housing_vulnerability_results_nn_match,
  any_housing_spell_between_18_19_subgroup_housing_vulnerability_results_full_match,
  any_housing_spell_between_18_19_subgroup_housing_vulnerability_results_nn_match,
  new_unsheltered_homelessness_between_18_19_subgroup_housing_vulnerability_results_full_match,
  new_unsheltered_homelessness_between_18_19_subgroup_housing_vulnerability_results_nn_match,
  new_ongoing_unsheltered_homelessness_between_18_19_subgroup_housing_vulnerability_results_full_match,
  new_ongoing_unsheltered_homelessness_between_18_19_subgroup_housing_vulnerability_results_nn_match,
  new_short_term_accommodation_required_between_18_19_subgroup_housing_vulnerability_results_full_match,
  new_short_term_accommodation_required_between_18_19_subgroup_housing_vulnerability_results_nn_match,
  new_ongoing_short_term_accommodation_required_between_18_19_subgroup_housing_vulnerability_results_full_match,
  new_ongoing_short_term_accommodation_required_between_18_19_subgroup_housing_vulnerability_results_nn_match,
  count_total_spells_18_19_subgroup_housing_vulnerability_results_full_match,
  count_total_spells_18_19_subgroup_housing_vulnerability_results_nn_match,
  time_in_housing_support_18_19_subgroup_housing_vulnerability_results_full_match,
  time_in_housing_support_18_19_subgroup_housing_vulnerability_results_nn_match
) |>
  rename(
    metrics = 1,
    in_housing_spell_on_date_18_prior_homelessness_full_match = in_housing_spell_on_date_18_prior_homelessness...2,                             
    in_housing_spell_on_date_18_no_prior_homelessness_full_match = in_housing_spell_on_date_18_no_prior_homelessness...3,                                
    in_housing_spell_on_date_18_prior_homelessness_nn_match = in_housing_spell_on_date_18_prior_homelessness...5,                                  
    in_housing_spell_on_date_18_no_prior_homelessness_nn_match = in_housing_spell_on_date_18_no_prior_homelessness...6,                                 
    in_housing_spell_on_date_19_prior_homelessness_full_match = in_housing_spell_on_date_19_prior_homelessness...8,                                 
    in_housing_spell_on_date_19_no_prior_homelessness_full_match = in_housing_spell_on_date_19_no_prior_homelessness...9,                                  
    in_housing_spell_on_date_19_prior_homelessness_nn_match = in_housing_spell_on_date_19_prior_homelessness...11,                                 
    in_housing_spell_on_date_19_no_prior_homelessness_nn_match = in_housing_spell_on_date_19_no_prior_homelessness...12,                                 
    new_housing_spell_between_18_19_prior_homelessness_full_match = new_housing_spell_between_18_19_prior_homelessness...14,                      
    new_housing_spell_between_18_19_no_prior_homelessness_full_match = new_housing_spell_between_18_19_no_prior_homelessness...15,                             
    new_housing_spell_between_18_19_prior_homelessness_nn_match = new_housing_spell_between_18_19_prior_homelessness...17,                       
    new_housing_spell_between_18_19_no_prior_homelessness_nn_match = new_housing_spell_between_18_19_no_prior_homelessness...18,                             
    any_housing_spell_between_18_19_prior_homelessness_full_match = any_housing_spell_between_18_19_prior_homelessness...20,                       
    any_housing_spell_between_18_19_no_prior_homelessness_full_match = any_housing_spell_between_18_19_no_prior_homelessness...21,                             
    any_housing_spell_between_18_19_prior_homelessness_nn_match = any_housing_spell_between_18_19_prior_homelessness...23,                       
    any_housing_spell_between_18_19_no_prior_homelessness_nn_match = any_housing_spell_between_18_19_no_prior_homelessness...24,                         
    new_unsheltered_homelessness_between_18_19_prior_homelessness_full_match = new_unsheltered_homelessness_between_18_19_prior_homelessness...26,            
    new_unsheltered_homelessness_between_18_19_no_prior_homelessness_full_match = new_unsheltered_homelessness_between_18_19_no_prior_homelessness...27,                  
    new_unsheltered_homelessness_between_18_19_prior_homelessness_nn_match = new_unsheltered_homelessness_between_18_19_prior_homelessness...29,            
    new_unsheltered_homelessness_between_18_19_no_prior_homelessness_nn_match = new_unsheltered_homelessness_between_18_19_no_prior_homelessness...30,                 
    new_ongoing_unsheltered_homelessness_between_18_19_prior_homelessness_full_match = new_ongoing_unsheltered_homelessness_between_18_19_prior_homelessness...32,    
    new_ongoing_unsheltered_homelessness_between_18_19_no_prior_homelessness_full_match = new_ongoing_unsheltered_homelessness_between_18_19_no_prior_homelessness...33,          
    new_ongoing_unsheltered_homelessness_between_18_19_prior_homelessness_nn_match = new_ongoing_unsheltered_homelessness_between_18_19_prior_homelessness...35,    
    new_ongoing_unsheltered_homelessness_between_18_19_no_prior_homelessness_nn_match = new_ongoing_unsheltered_homelessness_between_18_19_no_prior_homelessness...36,          
    new_short_term_accommodation_required_between_18_19_prior_homelessness_full_match = new_short_term_accommodation_required_between_18_19_prior_homelessness...38,   
    new_short_term_accommodation_required_between_18_19_no_prior_homelessness_full_match = new_short_term_accommodation_required_between_18_19_no_prior_homelessness...39,         
    new_short_term_accommodation_required_between_18_19_prior_homelessness_nn_match = new_short_term_accommodation_required_between_18_19_prior_homelessness...41,   
    new_short_term_accommodation_required_between_18_19_no_prior_homelessness_nn_match = new_short_term_accommodation_required_between_18_19_no_prior_homelessness...42,        
    new_ongoing_short_term_accommodation_required_between_18_19_prior_homelessness_full_match = new_ongoing_short_term_accommodation_required_between_18_19_prior_homelessness...44,
    new_ongoing_short_term_accommodation_required_between_18_19_no_prior_homelessness_full_match = new_ongoing_short_term_accommodation_required_between_18_19_no_prior_homelessness...45, 
    new_ongoing_short_term_accommodation_required_between_18_19_prior_homelessness_nn_match = new_ongoing_short_term_accommodation_required_between_18_19_prior_homelessness...47,
    new_ongoing_short_term_accommodation_required_between_18_19_no_prior_homelessness_nn_match = new_ongoing_short_term_accommodation_required_between_18_19_no_prior_homelessness...48, 
    count_total_spells_18_19_prior_homelessness_full_match = count_total_spells_18_19_prior_homelessness...50,                              
    count_total_spells_18_19_no_prior_homelessness_full_match = count_total_spells_18_19_no_prior_homelessness...51,                                    
    count_total_spells_18_19_prior_homelessness_nn_match = count_total_spells_18_19_prior_homelessness...53,                             
    count_total_spells_18_19_no_prior_homelessness_nn_match = count_total_spells_18_19_no_prior_homelessness...54,                                    
    time_in_housing_support_18_19_prior_homelessness_full_match = time_in_housing_support_18_19_prior_homelessness...56,                           
    time_in_housing_support_18_19_no_prior_homelessness_full_match = time_in_housing_support_18_19_no_prior_homelessness...57,                               
    time_in_housing_support_18_19_prior_homelessness_nn_match = time_in_housing_support_18_19_prior_homelessness...59,                          
    time_in_housing_support_18_19_no_prior_homelessness_nn_match = time_in_housing_support_18_19_no_prior_homelessness...60 
  ) |>
  select(
    -starts_with("reported_metrics...")
  )
    
# export results
write_csv(
  combined_subgroup_housing_vulnerability_results,
  "U:/pyi-paper/output/treatment_effect_results/combined_subgroup_housing_vulnerability_results.csv"
)

#-------------------------------------------------------------------------------
# 6. Table 6: Sensitivity analysis results
#-------------------------------------------------------------------------------

# read in all sensitivity results
sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/sensitivity_analysis/",
  common_suffix = "tip_results_full_match"
)

sherlock_reads_files_ending(
  folder_path = "U:/pyi-paper/output/sensitivity_analysis/",
  common_suffix = "tip_results_nn_match"
)

in_housing_spell_on_date_18_tip_results_full_match <- in_housing_spell_on_date_18_tip_results_full_match |> 
  mutate(
    outcome = "In homelessness spell on 18th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

in_housing_spell_on_date_18_tip_results_nn_match <- in_housing_spell_on_date_18_tip_results_nn_match |> 
  mutate(
    outcome = "In homelessness spell on 18th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

in_housing_spell_on_date_19_tip_results_full_match <- in_housing_spell_on_date_19_tip_results_full_match |> 
  mutate(
    outcome = "In homelessness spell on 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

in_housing_spell_on_date_19_tip_results_nn_match <- in_housing_spell_on_date_19_tip_results_nn_match |> 
  mutate(
    outcome = "In homelessness spell on 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

any_housing_spell_between_18_19_tip_results_full_match <- any_housing_spell_between_18_19_tip_results_full_match |> 
  mutate(
    outcome = "Any homelessness spell between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

any_housing_spell_between_18_19_tip_results_nn_match <- any_housing_spell_between_18_19_tip_results_nn_match |> 
  mutate(
    outcome = "Any homelessness spell between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_housing_spell_between_18_19_tip_results_full_match <- new_housing_spell_between_18_19_tip_results_full_match |> 
  mutate(
    outcome = "New homelessness spell between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

new_housing_spell_between_18_19_tip_results_nn_match <- new_housing_spell_between_18_19_tip_results_nn_match |> 
  mutate(
    outcome = "New homelessness spell between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_unsheltered_homelessness_between_18_19_tip_results_full_match <- new_unsheltered_homelessness_between_18_19_tip_results_full_match |> 
  mutate(
    outcome = "New unsheltered homelessness spell between 18th and 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

new_unsheltered_homelessness_between_18_19_tip_results_nn_match <- new_unsheltered_homelessness_between_18_19_tip_results_nn_match |> 
  mutate(
    outcome = "New unsheltered homelessness spell between 18th and 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_ongoing_unsheltered_homelessness_between_18_19_tip_results_full_match <- new_ongoing_unsheltered_homelessness_between_18_19_tip_results_full_match |> 
  mutate(
    outcome = "New or ongoing unsheltered homelessness spell between 18th and 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

new_ongoing_unsheltered_homelessness_between_18_19_tip_results_nn_match <- new_ongoing_unsheltered_homelessness_between_18_19_tip_results_nn_match |> 
  mutate(
    outcome = "New or ongoing unsheltered homelessness spell between 18th and 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_ongoing_short_term_accommodation_required_between_18_19_tip_results_full_match <- new_ongoing_short_term_accommodation_required_between_18_19_tip_results_full_match |> 
  mutate(
    outcome = "In new or ongoing homelessness spell that requires short term accommodation between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

new_ongoing_short_term_accommodation_required_between_18_19_tip_results_nn_match <- new_ongoing_short_term_accommodation_required_between_18_19_tip_results_nn_match |> 
  mutate(
    outcome = "In new or ongoing homelessness spell that requires short term accommodation between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_short_term_accommodation_required_between_18_19_tip_results_full_match <- new_short_term_accommodation_required_between_18_19_tip_results_full_match |> 
  mutate(
    outcome = "In new homelessness spell that requires short term accommodation between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Full"
  )

new_short_term_accommodation_required_between_18_19_tip_results_nn_match <- new_short_term_accommodation_required_between_18_19_tip_results_nn_match |> 
  mutate(
    outcome = "In new homelessness spell that requires short term accommodation between 18th & 19th birthday"
  ) |>
  rename(
    confounder_outcome_relationship = 4
  ) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

# combine results
combined_tip_analysis_results <- bind_rows(
    in_housing_spell_on_date_18_tip_results_full_match,
    in_housing_spell_on_date_18_tip_results_nn_match,
    in_housing_spell_on_date_19_tip_results_full_match,
    in_housing_spell_on_date_19_tip_results_nn_match,
    new_housing_spell_between_18_19_tip_results_full_match,
    new_housing_spell_between_18_19_tip_results_nn_match,
    any_housing_spell_between_18_19_tip_results_full_match,
    any_housing_spell_between_18_19_tip_results_nn_match,
    new_unsheltered_homelessness_between_18_19_tip_results_full_match,
    new_unsheltered_homelessness_between_18_19_tip_results_nn_match,
    new_ongoing_unsheltered_homelessness_between_18_19_tip_results_full_match,
    new_ongoing_unsheltered_homelessness_between_18_19_tip_results_nn_match,
    new_short_term_accommodation_required_between_18_19_tip_results_full_match,
    new_short_term_accommodation_required_between_18_19_tip_results_nn_match,
    new_ongoing_short_term_accommodation_required_between_18_19_tip_results_full_match,
    new_ongoing_short_term_accommodation_required_between_18_19_tip_results_nn_match
  ) |>
  select(
    outcome,
    matching_specification,
    scenario,
    confounder_prevalence_exposed,
    confounder_prevalence_unexposed,
    confounder_outcome_relationship 
  )

# export results
write_csv(
  combined_tip_analysis_results,
  "U:/pyi-paper/output/sensitivity_analysis/combined_tip_analysis_results.csv"
)

#-------------------------------------------------------------------------------
# 7. Figure 6: SMD plot data
#-------------------------------------------------------------------------------

# read modelling data
in_housing_spell_on_date_18_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/in_housing_spell_on_date_18_smd_results_full_match.RDS"
in_housing_spell_on_date_18_smd_results_full_match <- readRDS(in_housing_spell_on_date_18_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

in_housing_spell_on_date_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/in_housing_spell_on_date_19_smd_results_full_match.RDS"
in_housing_spell_on_date_19_smd_results_full_match <- readRDS(in_housing_spell_on_date_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

new_housing_spell_between_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/new_housing_spell_between_18_19_smd_results_full_match.RDS"
new_housing_spell_between_18_19_smd_results_full_match <- readRDS(new_housing_spell_between_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

any_housing_spell_between_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/any_housing_spell_between_18_19_smd_results_full_match.RDS"
any_housing_spell_between_18_19_smd_results_full_match <- readRDS(any_housing_spell_between_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

time_in_housing_support_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/time_in_housing_support_18_19_smd_results_full_match.RDS"
time_in_housing_support_18_19_smd_results_full_match <- readRDS(time_in_housing_support_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

new_unsheltered_homelessness_between_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/new_unsheltered_homelessness_between_18_19_smd_results_full_match.RDS"
new_unsheltered_homelessness_between_18_19_smd_results_full_match <- readRDS(new_unsheltered_homelessness_between_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

new_ongoing_unsheltered_homelessness_between_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/new_ongoing_unsheltered_homelessness_between_18_19_smd_results_full_match.RDS"
new_ongoing_unsheltered_homelessness_between_18_19_smd_results_full_match <- readRDS(new_ongoing_unsheltered_homelessness_between_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

new_short_term_accommodation_required_between_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/new_short_term_accommodation_required_between_18_19_smd_results_full_match.RDS"
new_short_term_accommodation_required_between_18_19_smd_results_full_match <- readRDS(new_short_term_accommodation_required_between_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

new_ongoing_short_term_accommodation_required_between_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/new_ongoing_short_term_accommodation_required_between_18_19_smd_results_full_match.RDS"
new_ongoing_short_term_accommodation_required_between_18_19_smd_results_full_match <- readRDS(new_ongoing_short_term_accommodation_required_between_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

count_total_spells_18_19_smd_results_full_match_file_location <- "U:/pyi-paper/output/plot_data/count_total_spells_18_19_smd_results_full_match.RDS"
count_total_spells_18_19_smd_results_full_match <- readRDS(count_total_spells_18_19_smd_results_full_match_file_location) |>
  mutate(
    matching_specification = "Full"
  )

in_housing_spell_on_date_18_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/in_housing_spell_on_date_18_smd_results_nn_match.RDS"
in_housing_spell_on_date_18_smd_results_nn_match <- readRDS(in_housing_spell_on_date_18_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

in_housing_spell_on_date_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/in_housing_spell_on_date_19_smd_results_nn_match.RDS"
in_housing_spell_on_date_19_smd_results_nn_match <- readRDS(in_housing_spell_on_date_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_housing_spell_between_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/new_housing_spell_between_18_19_smd_results_nn_match.RDS"
new_housing_spell_between_18_19_smd_results_nn_match <- readRDS(new_housing_spell_between_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

any_housing_spell_between_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/any_housing_spell_between_18_19_smd_results_nn_match.RDS"
any_housing_spell_between_18_19_smd_results_nn_match <- readRDS(any_housing_spell_between_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

time_in_housing_support_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/time_in_housing_support_18_19_smd_results_nn_match.RDS"
time_in_housing_support_18_19_smd_results_nn_match <- readRDS(time_in_housing_support_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_unsheltered_homelessness_between_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/new_unsheltered_homelessness_between_18_19_smd_results_nn_match.RDS"
new_unsheltered_homelessness_between_18_19_smd_results_nn_match <- readRDS(new_unsheltered_homelessness_between_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_ongoing_unsheltered_homelessness_between_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/new_ongoing_unsheltered_homelessness_between_18_19_smd_results_nn_match.RDS"
new_ongoing_unsheltered_homelessness_between_18_19_smd_results_nn_match <- readRDS(new_ongoing_unsheltered_homelessness_between_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_short_term_accommodation_required_between_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/new_short_term_accommodation_required_between_18_19_smd_results_nn_match.RDS"
new_short_term_accommodation_required_between_18_19_smd_results_nn_match <- readRDS(new_short_term_accommodation_required_between_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

new_ongoing_short_term_accommodation_required_between_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/new_ongoing_short_term_accommodation_required_between_18_19_smd_results_nn_match.RDS"
new_ongoing_short_term_accommodation_required_between_18_19_smd_results_nn_match <- readRDS(new_ongoing_short_term_accommodation_required_between_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

count_total_spells_18_19_smd_results_nn_match_file_location <- "U:/pyi-paper/output/plot_data/count_total_spells_18_19_smd_results_nn_match.RDS"
count_total_spells_18_19_smd_results_nn_match <- readRDS(count_total_spells_18_19_smd_results_nn_match_file_location) |>
  mutate(
    matching_specification = "Nearest Neighbour"
  )

# merge data
smd_plot_data <- bind_rows(
  in_housing_spell_on_date_18_smd_results_full_match,
  in_housing_spell_on_date_18_smd_results_nn_match,
  in_housing_spell_on_date_19_smd_results_full_match,
  in_housing_spell_on_date_19_smd_results_nn_match,
  new_housing_spell_between_18_19_smd_results_full_match,
  new_housing_spell_between_18_19_smd_results_nn_match,
  any_housing_spell_between_18_19_smd_results_full_match,
  any_housing_spell_between_18_19_smd_results_nn_match,
  new_unsheltered_homelessness_between_18_19_smd_results_full_match,
  new_unsheltered_homelessness_between_18_19_smd_results_nn_match,
  new_ongoing_unsheltered_homelessness_between_18_19_smd_results_full_match,
  new_ongoing_unsheltered_homelessness_between_18_19_smd_results_nn_match,
  new_short_term_accommodation_required_between_18_19_smd_results_full_match,
  new_short_term_accommodation_required_between_18_19_smd_results_nn_match,
  new_ongoing_short_term_accommodation_required_between_18_19_smd_results_full_match,
  new_ongoing_short_term_accommodation_required_between_18_19_smd_results_nn_match,
  time_in_housing_support_18_19_smd_results_full_match,
  time_in_housing_support_18_19_smd_results_nn_match,
  count_total_spells_18_19_smd_results_full_match,
  count_total_spells_18_19_smd_results_nn_match
) |>
  mutate(
    matching_specification = factor(matching_specification),
    outcome = factor(outcome),
    group = factor(group),
    strata = factor(strata)
  ) |>
  mutate(
    group = factor(
      group,
      levels = c(
        "Overall",
        "Sex",
        "Aboriginal",
        "Housing vulnerability")
    ),
    strata = factor(
      strata,
      levels = c(
        "Overall",
        "Male",
        "Female",
        "Aboriginal",
        "Non-Aboriginal",
        "Prior homelessness",
        "No prior homelessness"),
      ordered = TRUE
    )
  )

saveRDS(
  smd_plot_data,
  "U:/pyi-paper/output/plot_data/smd_plot_data.RDS"
)