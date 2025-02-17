#==============================================================================#
# RUN ALL MODELS                                                               #
#==============================================================================#

# outcome #1: in_housing_spell_on_date_18
source("./code/analysis/analysis_outcome_model_in_housing_spell_on_date_18_nn_match.R")
source("./code/analysis/analysis_outcome_model_in_housing_spell_on_date_18_full_match.R")

# outcome #2: in_housing_spell_on_date_19  
source("./code/analysis/analysis_outcome_model_in_housing_spell_on_date_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_in_housing_spell_on_date_19_full_match.R")

# outcome #3: new_housing_spell_between_18_19 
source("./code/analysis/analysis_outcome_model_new_housing_spell_between_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_new_housing_spell_between_18_19_full_match.R")

# outcome #4: any_housing_spell_between_18_19 
source("./code/analysis/analysis_outcome_model_any_housing_spell_between_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_any_housing_spell_between_18_19_full_match.R")

# outcome #5: time_in_housing_support_18_19  
source("./code/analysis/analysis_outcome_model_time_in_housing_support_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_time_in_housing_support_18_19_full_match.R")

# outcome #6: new_unsheltered_homelessness_between_18_19  
source("./code/analysis/analysis_outcome_model_new_unsheltered_homelessness_between_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_new_unsheltered_homelessness_between_18_19_full_match.R")

# outcome #7: new_ongoing_unsheltered_homelessness_between_18_19 
source("./code/analysis/analysis_outcome_model_new_ongoing_unsheltered_homelessness_between_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_new_ongoing_unsheltered_homelessness_between_18_19_full_match.R")

# outcome #8: new_short_term_accommodation_required_between_18_19   
source("./code/analysis/analysis_outcome_model_new_short_term_accommodation_required_between_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_new_short_term_accommodation_required_between_18_19_full_match.R")

# outcome #9: new_ongoing_short_term_accommodation_required_between_18_19  
source("./code/analysis/analysis_outcome_model_new_ongoing_short_term_accommodation_required_between_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_new_ongoing_short_term_accommodation_required_between_18_19_full_match.R")

# outcome #10: count_total_spells_18_19
source("./code/analysis/analysis_outcome_model_count_total_spells_18_19_nn_match.R")
source("./code/analysis/analysis_outcome_model_count_total_spells_18_19_full_match.R")

# visualisation data for export
source("./code/visualisation/visualisation_data_cleaning.R")
