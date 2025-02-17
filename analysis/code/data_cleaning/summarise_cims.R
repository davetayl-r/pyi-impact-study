# ========================================================================= #
# PYI Impact Study                                                          #
# Clean and summarise CIMS file                                             #
# Author: David Taylor                                                      #
# Date: 15/12/2024                                                          #
# ========================================================================= #

# This code clean and summarise cims file to: 
# a) find kids in pyi matching pool who are in cims, 
# b) summarise their interactions with shs
# note: file checks_summarise_cims_file.R performs checks to investigate the 
# integrity of this data and other assumptions

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages 
library(tidyverse)
library(readxl)
library(stringr)
library(lubridate)

# set up workspace preferences using custom helper function
workspace_setup()

# read oohc data for ids and dob
summarised_oohc_file_location <- "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_fifteen.RDS"
summarised_oohc_file <- readRDS(summarised_oohc_file_location)

# read cims file
cims_file_location <- "P:/psp/PSP/flat_files_final/CIMS/SHS_CIMS_all_MASTER.rds"
raw_cims_file <- readRDS(cims_file_location)

# get ids and dob from matching pool
matching_pool_ids <- summarised_oohc_file |>
  select(
    childstory_id,
    aihw_slk,
    date_of_birth,
    date_16,
    date_18,
    date_19
  ) |>
  ungroup()

# join cims data to matching pool
cims_data_summarise_step_one <- dplyr::inner_join( # perhaps change this to left_join post summarising
  matching_pool_ids,
  raw_cims_file,
    by = c("aihw_slk"),
  relationship = "many-to-many"
  ) |>
  # rename dob
  rename(
    date_of_birth = date_of_birth.x
  )

# identify if presentation reason is housing related
cims_data_summarise_step_two <- cims_data_summarise_step_one |>
  mutate(
  # ensure variables are numeric
    across(
      c(
        residential_dwelling_wkbef,
        residential_dwelling_present,
        Previously_Homeless_Mth1, 
        Previously_Homeless_Mth2,
        Short_Term_Accom_ind),
      ~as.numeric(replace_na(.,0))),
    # primary homelessness indicators
    unsheltered_homelessness_current = case_when(
      # current residential dwelling is a tent (=1), improvised dwelling (=6), no dwelling (in the open)(=7) or motor vehicle (=8)
      residential_dwelling_present %in% c(1, 6, 7, 8) ~ 1,
      TRUE ~ 0
      ),
    unsheltered_homelessness_last_week = case_when(
      # residential dwelling in last week was a tent (=1), improvised dwelling (=6), no dwelling (in the open)(=7) or motor vehicle (=8)
      residential_dwelling_wkbef %in% c(1, 6, 7, 8) |
      # sleeping rough or in non-conventional accommodation in last week
      Previously_Homeless_Mth1 == 1 ~ 1,
      TRUE ~ 0
      ),
    # service need indicators
    requires_housing_assistance = case_when( 
      # residential dwelling in last week was a tent (=1), improvised dwelling (=6), no dwelling (in the open)(=7) or motor vehicle (=8)
      residential_dwelling_wkbef %in% c(1, 6, 7, 8) |
        # current residential dwelling is a tent (=1), improvised dwelling (=6), no dwelling (in the open)(=7) or motor vehicle (=8)
        residential_dwelling_present %in% c(1, 6, 7, 8) |
        # sleeping rough or in non-conventional accommodation in last week
        Previously_Homeless_Mth1 == 1 |
        # in short-term or emergency accommodation in the last week
        Previously_Homeless_Mth2 == 1 |
        # require short-term accommodation 
        Short_Term_Accom_ind == 1 ~ 1,
      TRUE ~ 0
      ),
    # use same criteria as original analysis
    requires_housing_assistance_original = case_when( 
      # sleeping rough or in non-conventional accommodation in last week
      Previously_Homeless_Mth1 == 1 |
      # in short-term or emergency accommodation in the last week
      Previously_Homeless_Mth2 == 1 |
      # require short-term accommodation 
      Short_Term_Accom_ind == 1 ~ 1,
      TRUE ~ 0
      ),
    # identify presentations for non housing reasons
    non_housing_presentations = case_when(
      unsheltered_homelessness_current == 0 &
      unsheltered_homelessness_last_week == 0 &
      requires_housing_assistance == 0 ~ 1,
      TRUE ~ 0
    )
  ) |>
  # drop cases where presentation reason is not housing-related
  filter(
    non_housing_presentations == 0
  )

# identify spells
cims_data_summarise_step_three <- cims_data_summarise_step_two |>
  select(
    aihw_slk,
    childstory_id,
    support_period_start_date,
    support_period_end_date,
    date_of_birth,
    date_16,
    date_18,
    date_19,
    unsheltered_homelessness_current,
    unsheltered_homelessness_last_week,
    requires_housing_assistance
  ) |>
  # ensure variables are dates
  mutate(
    support_period_start_date = as.Date(support_period_start_date),
    support_period_end_date = as.Date(support_period_end_date)
  ) |>
  # arrange data
  arrange(
    aihw_slk,
    support_period_start_date,
    desc(support_period_end_date)
  ) |>
  # add columns to store information
  mutate(
    sp_type = as.character(NA),
    spell_start_date = as.Date(NA),
    spell_end_date = as.Date(NA),
    spell_number = as.numeric(NA),
    spell_id = as.numeric(NA)
  ) 

# initialise counters
sequence <- 0
sequence_number <- 0

# run for loop to initialise first record
for (i in 1:1) {
  cims_data_summarise_step_three$sp_type[i] <- "new_spell" 
  cims_data_summarise_step_three$spell_start_date[i] <- cims_data_summarise_step_three$support_period_start_date[i]
  cims_data_summarise_step_three$spell_end_date[i] <- cims_data_summarise_step_three$support_period_end_date[i] 
  sequence_number <- sequence_number + 1
  cims_data_summarise_step_three$spell_number[i] <- sequence_number  
  sequence <- sequence + 1
  cims_data_summarise_step_three$spell_id[i] <- sequence
}

for (i in 2:nrow(cims_data_summarise_step_three)){
  if(
    cims_data_summarise_step_three$aihw_slk[i]!=cims_data_summarise_step_three$aihw_slk[i-1]){
    # find first appearance for individual
    cims_data_summarise_step_three$sp_type[i] <- "new_spell" 
    cims_data_summarise_step_three$spell_start_date[i] <- cims_data_summarise_step_three$support_period_start_date[i]
    cims_data_summarise_step_three$spell_end_date[i]<- cims_data_summarise_step_three$support_period_end_date[i] 
    sequence_number <- 1
    cims_data_summarise_step_three$spell_number[i] <- sequence_number
    sequence <- sequence + 1
    cims_data_summarise_step_three$spell_id[i] <- sequence
  } 
  else if (
    cims_data_summarise_step_three$aihw_slk[i]==cims_data_summarise_step_three$aihw_slk[i-1] & 
    cims_data_summarise_step_three$support_period_start_date[i]  <= cims_data_summarise_step_three$spell_end_date[i-1] & 
    cims_data_summarise_step_three$support_period_start_date[i] >= cims_data_summarise_step_three$spell_start_date[i-1] & 
    cims_data_summarise_step_three$support_period_end_date[i] <= cims_data_summarise_step_three$spell_end_date[i-1]
  ){
    # find out if the same record for that individual is contained within the previous record
    cims_data_summarise_step_three$sp_type[i] <- "contained_within" 
    cims_data_summarise_step_three$spell_start_date[i] <- cims_data_summarise_step_three$spell_start_date[i-1]
    cims_data_summarise_step_three$spell_end_date[i] <- cims_data_summarise_step_three$spell_end_date[i-1] 
    cims_data_summarise_step_three$spell_number[i] <- sequence_number
    cims_data_summarise_step_three$spell_id[i] <- sequence
  } 
  else if (
    cims_data_summarise_step_three$aihw_slk[i] == cims_data_summarise_step_three$aihw_slk[i-1] & 
    cims_data_summarise_step_three$support_period_start_date[i]  <= cims_data_summarise_step_three$spell_end_date[i-1] & 
    cims_data_summarise_step_three$support_period_start_date[i] >= cims_data_summarise_step_three$spell_start_date[i-1] & 
    cims_data_summarise_step_three$support_period_end_date[i] > cims_data_summarise_step_three$spell_end_date[i-1])  
  {
    # find out if record is a transfer
    cims_data_summarise_step_three$sp_type[i] <- "transferred_to" 
    cims_data_summarise_step_three$spell_start_date[i] <- cims_data_summarise_step_three$spell_start_date[i-1]
    cims_data_summarise_step_three$spell_number[i] <- sequence_number
    cims_data_summarise_step_three$spell_id[i] <- sequence
    cims_data_summarise_step_three$spell_end_date[i] <- cims_data_summarise_step_three$support_period_end_date[i]
  } 
  else if (
    cims_data_summarise_step_three$aihw_slk[i] == cims_data_summarise_step_three$aihw_slk[i-1] & 
    cims_data_summarise_step_three$support_period_start_date[i] > cims_data_summarise_step_three$spell_end_date[i-1]) 
  {
    # find out if record is new
    cims_data_summarise_step_three$sp_type[i] <- "new_spell" 
    cims_data_summarise_step_three$spell_start_date[i] <- cims_data_summarise_step_three$support_period_start_date[i]
    cims_data_summarise_step_three$spell_end_date[i]<- cims_data_summarise_step_three$support_period_end_date[i]
    sequence_number <- sequence_number + 1
    cims_data_summarise_step_three$spell_number[i] <- sequence_number 
    sequence <- sequence + 1
    cims_data_summarise_step_three$spell_id[i] <- sequence
  } 
  else {cims_data_summarise_step_three$sp_type[i] <- "error"}
}

# summarise data per spell 
cims_data_summarise_step_four <- cims_data_summarise_step_three |>   
  # group by individual and spell
  group_by(
    aihw_slk,
    childstory_id,
    date_of_birth,
    date_16,
    date_18,
    date_19,
    spell_id,
    spell_start_date
  ) |>
  # summarise data by spell
  summarise(
    unsheltered_homelessness_current = max(unsheltered_homelessness_current),
    unsheltered_homelessness_last_week = max(unsheltered_homelessness_last_week),
    requires_housing_assistance = max(requires_housing_assistance),
    spell_end_date = max(spell_end_date)
  ) |>
  # select vars
  ungroup() |>
  select(
    aihw_slk,
    childstory_id,
    date_of_birth,
    date_16,
    date_18,
    date_19,
    spell_id,
    spell_start_date,
    spell_end_date,
    unsheltered_homelessness_current,
    unsheltered_homelessness_last_week,
    requires_housing_assistance
  )

# create variables
cims_data_summarise_step_five <- cims_data_summarise_step_four |>
  # group by individual
  group_by(
    aihw_slk
  ) |>
  # create variables for outcome measurement and balance checks
  mutate(
    # any type of housing spell before 18th birthday
    any_housing_spell_before_18 = case_when(
      spell_start_date < date_18 ~ 1,
      TRUE ~ 0),
    # in housing spell on date_18
    in_housing_spell_on_date_18 = case_when(
      date_18 %within%
        interval(
          spell_start_date,
          spell_end_date
        ) ~ 1,
      TRUE ~ 0
    ),
    # in housing spell on date_19
    in_housing_spell_on_date_19 = case_when(
      date_19 %within%
        interval(
          spell_start_date,
          spell_end_date
        ) ~ 1,
      TRUE ~ 0
    ),
    # new housing spell between 18th and 19th birthday
    new_housing_spell_between_18_19 = case_when(
      spell_start_date %within% 
        interval(
          date_18,
          date_19
        ) ~ 1,
      TRUE ~ 0
    ),
    # any housing spell between 18th and 19th birthday
    any_housing_spell_between_18_19 = case_when(
      new_housing_spell_between_18_19 == 1 |
      in_housing_spell_on_date_18 == 1 ~ 1,
      TRUE ~ 0
    ),
    # calculate length of time in housing support between date_18 and date_19
    spell_start_date_18_19_censored = case_when(
      in_housing_spell_on_date_18 == 1 ~ date_18,
      TRUE ~ spell_start_date
      ),
    spell_end_date_18_19_censored = case_when(
      (new_housing_spell_between_18_19 == 1 | in_housing_spell_on_date_18 == 1) & spell_end_date > date_19 ~ date_19,
      TRUE ~ spell_end_date
      ), 
    time_in_housing_support_18_19 = case_when(
      (new_housing_spell_between_18_19 == 1 | in_housing_spell_on_date_18 == 1) ~ interval(spell_start_date_18_19_censored, spell_end_date_18_19_censored) %/% days(),
      TRUE ~ 0
      ),
    # new spell where individual reported being in unsheltered homelessness between 18th and 19th birthday
    new_unsheltered_homelessness_between_18_19 = case_when(
      new_housing_spell_between_18_19 == 1 & 
        (unsheltered_homelessness_current == 1 |
           unsheltered_homelessness_last_week) ~ 1,
      TRUE ~ 0
      ),
    # new or ongoing spell where individual reported being in unsheltered homelessness on or between 18th and 19th birthday
    new_ongoing_unsheltered_homelessness_between_18_19 = case_when(
      new_unsheltered_homelessness_between_18_19 == 1 |
        in_housing_spell_on_date_18 == 1 & 
        (unsheltered_homelessness_current == 1 |
           unsheltered_homelessness_last_week) ~ 1,
      TRUE ~ 0),
    # new spell where individual requires short term accommodation between 18th and 19th birthday
    new_short_term_accommodation_required_between_18_19 = case_when(
      new_housing_spell_between_18_19 == 1 & 
        requires_housing_assistance == 1 ~ 1,
      TRUE ~ 0),
    # new or ongoing spell where individual requires short term accommodation on or between 18th and 19th birthday
    new_ongoing_short_term_accommodation_required_between_18_19 = case_when(
      new_short_term_accommodation_required_between_18_19 == 1 |
        in_housing_spell_on_date_18 == 1 & 
        (requires_housing_assistance == 1) ~ 1,
      TRUE ~ 0),
    # new spell where individual requires short term accommodation between 18th and 19th birthday
    new_requires_housing_assistance_between_18_19 = case_when(
      new_housing_spell_between_18_19 == 1 & 
      requires_housing_assistance == 1 ~ 1,
      TRUE ~ 0),
    # new or ongoing spell where individual requires short term accommodation on or between 18th and 19th birthday
    new_ongoing_requires_housing_assistance_between_18_19 = case_when(
      new_requires_housing_assistance_between_18_19 == 1 |
      in_housing_spell_on_date_18 == 1 & 
      (requires_housing_assistance == 1) ~ 1,
      TRUE ~ 0),
    # count number of spells between date_18 and date_19
    count_total_spells_18_19 = sum(
      case_when(
        # in housing spell on date_18
        in_housing_spell_on_date_18 == 1 ~ 1,
        # new housing spell between 18th and 19th birthday
        new_housing_spell_between_18_19 == 1 ~ 1,
        TRUE ~ 0)
      )
    ) |>
  ungroup()

# summarise variables
cims_data_summarise_step_six <- cims_data_summarise_step_five |>
  group_by(
    aihw_slk,
    childstory_id
  ) |>
  summarise(
    any_housing_spell_before_18 = max(any_housing_spell_before_18),
    in_housing_spell_on_date_18 = max(in_housing_spell_on_date_18),
    in_housing_spell_on_date_19 = max(in_housing_spell_on_date_19),
    new_housing_spell_between_18_19 = max(new_housing_spell_between_18_19),
    any_housing_spell_between_18_19 = max(any_housing_spell_between_18_19),
    time_in_housing_support_18_19 = sum(time_in_housing_support_18_19),
    new_unsheltered_homelessness_between_18_19 = max(new_unsheltered_homelessness_between_18_19),
    new_ongoing_unsheltered_homelessness_between_18_19 = max(new_ongoing_unsheltered_homelessness_between_18_19),
    new_short_term_accommodation_required_between_18_19 = max(new_short_term_accommodation_required_between_18_19),
    new_ongoing_short_term_accommodation_required_between_18_19 = max(new_ongoing_short_term_accommodation_required_between_18_19),
    new_requires_housing_assistance_between_18_19 = max(new_requires_housing_assistance_between_18_19),
    new_ongoing_requires_housing_assistance_between_18_19 = max(new_ongoing_requires_housing_assistance_between_18_19), 
    count_total_spells_18_19 = max(count_total_spells_18_19)
  )

# export cims data
saveRDS(cims_data_summarise_step_six, "P:/pyi/pyi_update/data/processed_data/cims_data_summarise_step_six.RDS")
