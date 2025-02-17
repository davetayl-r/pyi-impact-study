# clean and summarise oohc file to calculate covariates for matching model

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages using custom helper function
required_packages <- c(
  "tidyverse",
  "readxl",
  "stringr",
  "lubridate",
  "progress"
)

snapshot_date <- "2024-11-01"

check_install_groundhog(
  required_packages,
  snapshot_date
)

# set up workspace preferences using custom helper function
workspace_setup()

# read rds files using custom function
sherlock_reads_files(
  folder_path = "P:/pyi/pyi_update/data/processed_data/",
  common_prefix = "oohc_data_summarise_step"
)

# check for conflicting records of gender
ids_inconsistent_gender_records <- check_record_consistency(
  oohc_data_summarise_step_one,
  "childstory_id",
  "Gender"
)

# A tibble: 0 x 3
# i 3 variables: childstory_id <chr>, count_distinct_categories <int>, categories <chr>

# check tally of gender, which can be recorded in one of five ways: 
# a) male, b) female, c) intersex or indeterminate, d) unknown, e) Not stated/inadequately described
gender_summary_oohc_file <- oohc_data_summarise_step_one |>
  select(
    childstory_id,
    Gender
  ) |>
  distinct() |>
  mutate(
    Gender = factor(Gender)
  ) |>
  group_by(
    Gender
  ) |>
  summarise(
    count = n()
  )

#Gender                            count
#<fct>                             <int>
#1 Female                            16789
#2 Intersex or indeterminate            17 - enough are within age range to warrant a category
#3 Male                              17740
#4 Not stated/inadequately described     4 - one is within age range, suggest combining with intersex or indeterminate category
#5 Unknown                               2 - we can ignore these - these are not in our age range

# check for conflicting records of aboriginality
ids_inconsistent_aboriginal_records <- check_record_consistency(
  oohc_data_summarise_step_one,
  "childstory_id",
  "AboriginalStatusFaCSGrouped"
)

# A tibble: 0 x 3
# i 3 variables: childstory_id <chr>, count_distinct_categories <int>, categories <chr>

# check tally of aboriginal, which can be recorded in one of four ways: 
# a) "", b) Aboriginal, c) Non-Aboriginal, d) Not Stated
aboriginal_summary_oohc_file <- oohc_data_summarise_step_one |>
  select(
    childstory_id,
    AboriginalStatusFaCSGrouped
  ) |>
  distinct() |>
  mutate(
    AboriginalStatusFaCSGrouped = factor(AboriginalStatusFaCSGrouped)
  ) |>
  group_by(
    AboriginalStatusFaCSGrouped
  ) |>
  summarise(
    count = n()
  )

#AboriginalStatusFaCSGrouped count
#<fct>                       <int>
#1 ""                              3 - one record within potential age range
#2 "Aboriginal"                13096
#3 "Non-Aboriginal"            21438
#4 "Not Stated"                   15 - >50% clearly out of potential age range
# missing data to be treated as "Non-Aboriginal"

# check to see if all kids who got pyi were under the personal responsibility of the minister
pyi_prm_test <- oohc_data_summarise_step_one |>
  filter(
    pyi_flag == 1
  ) |>
  mutate(
    not_oohc = case_when(
      PRGrouped %in% c("No Legal Orders", "Relative/Kinship Care - No Order", "Guardianship") ~ 1,
      TRUE ~ 0)
  ) |>
  mutate(
    PRM_any_time = case_when(
      PRGrouped %in% c("PR to Minister", "PR to Minister (Protected Person)") ~ 1,
      TRUE ~ 0)
  ) |>
  group_by(childstory_id) |>
  mutate(
    filter_var = max(not_oohc),
    prm_var = max(PRM_any_time)) |>
  filter(prm_var == 0) |>
  select(
    childstory_id,
    PRGrouped
  ) |>
  distinct(childstory_id)

# three individuals are not PRM during the evaluation period

# childstory_id
# <chr>        
# 1 C-00144446   
# 2 C-00181665   
# 3 C-01306538

# look for missing values in dates
sum(is.na(oohc_data_summarise_step_one$date_of_birth)) #9
sum(is.na(oohc_data_summarise_step_one$care_category_start_date)) #0
sum(is.na(oohc_data_summarise_step_one$care_category_end_date)) #0
sum(is.na(oohc_data_summarise_step_one$placement_start_date)) #0
sum(is.na(oohc_data_summarise_step_one$placement_end_date)) #0

# check all pyi kids are included once age eligible restrictions are applied
oohc_data_summarise_step_two |>
  filter(
    pyi_flag == 1
  ) |>
  select(
    childstory_id,
    date_of_birth,
    age_eligible_during_evaluation) |>
  distinct(
    childstory_id
  ) |>
  dim() # expecting 295

# check output from flagged_oohc_file

# visual check to see if check time in care flag works as expected
check_time_in_care_behaves <- oohc_data_summarise_step_six |>
  select(
    childstory_id,
    care_category_start_date,
    care_category_end_date,
    placement_start_date,
    placement_end_date_adjusted,
    time_in_care_start_eligible_window,
    time_in_care_end_eligible_window,
    eligible_time_in_care
  ) |>
  View()

# visual check to see if placement instability flag works as expected
check_placement_instabilty_behaves <- oohc_data_summarise_step_five |>
  select(
    childstory_id,
    care_category_start_date,
    care_category_end_date,
    placement_start_date,
    placement_end_date,
    date_benchmark,
    ActualPlacementType,
    duplicate_placement, 
    placement_more_thirty_days_not_duplicate, 
    placement_instability_in_oohc_placement, 
    placement_instability_in_last_six_months, 
    placement_instability_in_last_twelve_months, 
    count_placements_within_care_category, 
    count_placements_within_six_months, 
    count_placements_within_twelve_months, 
    placement_instability_binary 
  )

# visual check to see if permanent placement flag works as expected
check_permanent_placement_behaves <- oohc_data_summarise_step_five |>
  select(
    childstory_id,
    care_category_start_date,
    care_category_end_date,
    placement_start_date,
    placement_end_date,
    date_benchmark,
    PriorityPlacementPurpose,
    permanent_placement,
    permanent_placement_current
  ) |>
  filter(
    PriorityPlacementPurpose == "Permanent Care")

# visual check to see if resi care flag works as expected
check_residential_care_behaves <- oohc_data_summarise_step_five |>
  select(
    childstory_id,
    care_category_start_date,
    care_category_end_date,
    placement_start_date,
    placement_end_date,
    date_benchmark,
    ActualPlacementType,
    residential_care, 
    residential_care_current 
  ) |>
  filter(
    ActualPlacementType == "Residential Care")

# visual check to see if other vars work as expected
check_other_vars_behave <- oohc_data_summarise_step_five |>
  select(
    childstory_id,
    care_category_start_date,
    care_category_end_date,
    placement_start_date,
    placement_end_date,
    date_benchmark,
    PRGrouped,
    permanent_responsibility_minister, 
    permanent_responsibility_minister_current, 
    ActualPlacementType,
    self_placed,
    self_placed_current, 
    self_placed_after_16
  ) |>
  filter(
    ActualPlacementType %in% c(
      "Self Placed - Not Authorised",       
      "Self Placed - Not Authorised - Independent",
      "Self Placed - Not Authorised - Other Person/s",
      "Self Placed - Not Authorised - Parent/s",
      "Absent - Location Unknown") |
      PRGrouped == "PR to Minister"
  )

# find out how many kids who received pyi are tagged as eligible
oohc_data_summarise_step_eleven |>
  filter(
    pyi_flag == 1 &
    meets_pyi_eligibility_criteria_during_evaluation == 1
  ) |>
  distinct(
    childstory_id
  ) |>
  dim() #295

# find childstory_ids who are in step_eleven
ids_in_data <- step_eleven |>
  filter(
    pyi_flag == 1 &
    meets_pyi_eligibility_criteria_during_evaluation == 1
  ) |>
  distinct(
    childstory_id
  )

# find childstory_ids of individuals who got pyi but who are not in step_eleven
missing_ids <- oohc_data_summarise_step_one |>
  filter(
    pyi_flag == 1,
    !childstory_id %in% ids_in_data$childstory_id
  ) |>
  distinct(
    childstory_id
  ) |>
  dim() #5

# list missing ids
step_one |>
  filter(
    pyi_flag == 1,
    !childstory_id %in% ids_in_data$childstory_id
  ) |>
  select(
    childstory_id 
  ) |>
  distinct()

# inspect each individually
oohc_data_summarise_step_eleven |>
  filter(
    childstory_id == "C-00235726"
  ) |>
  View()


#1 C-00235726 M SWS PRM Placement ended age 16, ticks placement instability box but current placement appears to have ended without restoration  
#2 C-01683767 F WSW PRM Another placement closing before eligible without restoration, says 'planned move' in care   
#4 C-02581185 M WSW Young person incarcerated
#5 C-02325470 M NE PRM Placement recorded as ended due to 'Disruption Involving CYP', legal orders till open etc.
#6 C-02115156 M NBM PRM Restored to parents - should be ineligible

check_age_care_category_ended <- oohc_data_summarise_step_eleven |>
  filter(
    pyi_flag == 1,
    !childstory_id %in% ids_in_data$childstory_id
  ) |>
  group_by(
    childstory_id
  ) |>
  mutate(
    placement_count = row_number(),
    age_care_category_end_date = lubridate::time_length(
      interval(
        date_of_birth, 
        care_category_end_date), "years")
  ) |>
  filter(
    placement_count == max(placement_count)
  ) |>
  select(
    childstory_id,
    age_care_category_end_date
  )

# childstory_id age_care_category_end_date
# Groups:   childstory_id [8]
#childstory_id age_care_category_end_date
#<chr>                              <dbl>
#1 C-01501715                          15.2
#2 C-00784854                          16.1
#3 C-02325470                          16.5
#4 C-00235726                          15.7
#5 C-02115156                          15.9
#6 C-01683767                          15.9
#7 C-02581185                          15.6
#8 C-00528614                          17.9

# check how many folks "quality for pyi" under these conditions
oohc_data_summarise_step_eleven |>
  select(
    childstory_id,
    meets_pyi_eligibility_criteria_during_evaluation
  ) |>
  group_by(
    childstory_id
  ) |>
  mutate(
    meets_pyi_eligibility_criteria_during_evaluation = max(meets_pyi_eligibility_criteria_during_evaluation)
  ) |>
  distinct() |>
  dim() #4884

# how many kids are on this list
oohc_data_summarise_step_fifteen |>
  select(
    childstory_id
  ) |>
  distinct() |>
  dim() #1130

# check for duplicates
oohc_data_summarise_step_fifteen |>
  dplyr::count(
    childstory_id
  ) |>
  filter(
    n > 1
  ) # no duplicates 

# double check that all kids are on this list
oohc_data_summarise_step_fifteen |>
  filter(
    pyi_flag == 1
  ) |>
  select(
    childstory_id
  ) |>
  distinct() |>
  dim() #295

# check how many in the matching pool in step_one
oohc_data_summarise_step_one |>
  filter(
    matching_pool_flag == 1
  ) |>
  select(
    childstory_id
  ) |>
  distinct() |>
  dim() # 1299 

# check how many of them are in the updated data
oohc_data_summarise_step_fifteen |>
  filter(
    matching_pool_flag == 1
  ) |>
  select(
    childstory_id
  ) |>
  distinct() |>
  dim() # 483 of the original matching pool are here

483/1299 = 0.37%