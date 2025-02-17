# clean and summarise oohc file to calculate covariates for matching and outcome model
# note: file checks_summarise_oohc_file.R performs checks to investigate the integrity of this data and other assumptions

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages using custom helper function
required_packages <- c(
  "tidyverse",
  "readxl",
  "stringr",
  "lubridate"
)

snapshot_date <- "2024-11-01"

check_install_groundhog(
  required_packages,
  snapshot_date
)

# set up workspace preferences using custom helper function
workspace_setup()

# read data
oohc_file_location <- "P:/pyi/pyi_update/data/processed_data/oohc_file_pyi_flags.RDS"
oohc_file <- readRDS(oohc_file_location)

# set parameters for dates
#benchmark_age <- years(17) + months(6)
pyi_age_eligible_floor <- years(16) + months(9)
pyi_age_eligible_ceiling <- years(17) + months(6)
pyi_evaluation_start <- ymd("2018-01-01")
pyi_evaluation_end <- ymd("2021-01-01")

# step one: clean date variables and drop irrelevant columns
oohc_data_summarise_step_one <- oohc_file |>
  mutate(
    # convert dates from character vectors to dates
    date_of_birth = mdy(DOB),
    placement_start_date = ymd(PriorityPlacementStartDate),
    placement_end_date = ymd(PriorityPlacementEndDate),
    care_category_start_date = ymd(CareCategoryStartDate),
    care_category_end_date = ymd(CareCategoryEndDate),
    legal_orders_start_date = ymd(LegalOrderStartDate),
    legal_orders_end_date = ymd(LegalOrderEndDate),
    # create dates for range of events
    date_19 = add_with_rollback(date_of_birth, years(19)),
    date_18 = add_with_rollback(date_of_birth, years(18)),
    date_16 = add_with_rollback(date_of_birth, years(16))
  ) |>
  # create new end date variable for permanent placements with closure dates before age 18
  group_by(
    childstory_id
  ) |>
  arrange(
    placement_start_date
  ) |>
  mutate(
    # replace missing placement_end_dates with care_category_end_dates
    placement_end_date_adjusted = case_when(
      is.na(placement_end_date) ~ care_category_end_date,
      TRUE ~ placement_end_date)
  ) |>
  ungroup() |>
  # drop redundant columns
  select(
    -DOB,
    -PriorityPlacementStartDate,
    -PriorityPlacementEndDate,
    -ActualPlacementStartDate,
    -ActualPlacementEndDate,
    -CareCategoryStartDate,
    -CareCategoryEndDate,
    -LegalOrderStartDate,
    -LegalOrderEndDate,
    -ActualPlacementID,
    -Agency_ID,
    -PriorityPlacementID,
    -PriorityPlacementProviderGrouped,
    -PriorityPlacementAboriginalPlacementPrincipleFlag,
    -PRUngrouped
  )

saveRDS(oohc_data_summarise_step_one, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_one.RDS")
oohc_data_summarise_step_one <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_one.RDS")

# step two: drop individuals who are well outside age eligibility during evaluation period (mainly to speed up data processing)
oohc_data_summarise_step_two <- oohc_data_summarise_step_one |>
  group_by(
    childstory_id
  ) |>
  mutate(
    # create dates when individual enters and exits eligibility window for pyi
    eligible_age_window_start = add_with_rollback(
      date_of_birth,  pyi_age_eligible_floor),
    eligible_age_window_end = add_with_rollback(
      date_of_birth, pyi_age_eligible_ceiling),
    # flag if age-eligible during evaluation period
    age_eligible_during_evaluation = case_when(
      int_overlaps(
        interval(
          eligible_age_window_start,
          eligible_age_window_end),
        interval(
          pyi_evaluation_start,
          pyi_evaluation_end)
      ) ~ 1,
      TRUE ~ 0)
  ) |>
  filter(
    # drop cases with missing date_of_birth
    !is.na(date_of_birth),
    # keep cases where individual is age eligible during evaluation
    age_eligible_during_evaluation == 1,
    # drop cases where placement end date is before their DOB - there is one known case with a placement that opened and closed 10 years before their DOB
    !placement_end_date_adjusted < date_of_birth
    ) |>
  ungroup()

saveRDS(oohc_data_summarise_step_two, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_two.RDS")
oohc_data_summarise_step_two <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_two.RDS")

# step three: create flags for demographic details
oohc_data_summarise_step_three <- oohc_data_summarise_step_two |>
  mutate(
    # create flag for male (vs. not male)
    male = case_when(
      Gender == "Male" ~ 1,
      TRUE ~ 0),
    # create flag for female (vs. not female)
    female = case_when(
      Gender == "Female" ~ 1,
      TRUE ~ 0),
    # create flag for intersex, indeterminate or not stated (vs. not)
    intersex_other = case_when(
      Gender == "Intersex or indeterminate" |
        Gender == "Not stated/inadequately described" ~ 1,
      TRUE ~ 0),
    # create flag for aboriginal (vs. not)
    aboriginal = case_when(
      AboriginalStatusFaCSGrouped == "Aboriginal" ~ 1,
      TRUE ~ 0)
  ) |>
  # drop redundant demographic vars
  select(
    -Gender,
    -AboriginalStatusFaCSGrouped
  )

saveRDS(oohc_data_summarise_step_three, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_three.RDS")
oohc_data_summarise_step_three <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_three.RDS")

# step four: create flags for placement types within benchmark period
oohc_data_summarise_step_four <- oohc_data_summarise_step_three |>
  mutate(
    # create flag for residential care
    residential_care_placement = case_when(
      ActualPlacementType == "Residential Care" ~ 1,
      TRUE ~ 0),
    # create flag for residential care during eligibility window
    eligible_residential_care_placement = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          eligible_age_window_start,
          eligible_age_window_end)
        ) & residential_care_placement == 1 ~ 1,
      TRUE ~ 0), 
    # create flag for foster care placement
    foster_care_placement = case_when(
      ActualPlacementType == "Carer - Foster Carer" ~ 1,
      TRUE ~ 0),
    # create flag for foster care placement during eligibility window
    eligible_foster_care_placement = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          eligible_age_window_start,
          eligible_age_window_end)
      ) & foster_care_placement == 1 ~ 1,
      TRUE ~ 0),
    # create flag for kinship care placement
    kinship_care_placement = case_when(
      ActualPlacementType == "Carer - Relative or Kinship Carer" ~ 1,
      TRUE ~ 0),
    # create flag for kinship care placement during eligibility window
    eligible_kinship_care_placement = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          eligible_age_window_start,
          eligible_age_window_end)
      ) & kinship_care_placement == 1 ~ 1,
      TRUE ~ 0),
    # create flag for permanent placement 
    permanent_placement = case_when(
      PriorityPlacementPurpose == "Permanent Care" ~ 1,
      TRUE ~ 0),
    # create flag for current permanent placement during eligibility window
    eligible_permanent_placement = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          eligible_age_window_start,
          eligible_age_window_end)
      ) & permanent_placement == 1 ~ permanent_placement,
      TRUE ~ 0),
    # create flag for permanent responsibility to minister
    permanent_responsibility_minister = case_when(
      PRGrouped == "PR to Minister" ~ 1,
      PRGrouped == "PR to Minister (Protected Person)" ~ 1,
      TRUE ~ 0),
    # create flag for current permanent responsibility to minister during eligibility window
    permanent_responsibility_minister_eligiblity_period = case_when(
      int_overlaps(
        interval(
          legal_orders_start_date, 
          legal_orders_end_date), 
        interval(
          eligible_age_window_start,
          eligible_age_window_end)
      ) & permanent_responsibility_minister == 1 ~ 1,
      TRUE ~ 0),
    # create flag for self placed and/or runaway
    self_placed_missing_absent_from_placement = case_when(
      ActualPlacementType %in% c(
        "Self Placed - Not Authorised",       
        "Self Placed - Not Authorised - Independent",
        "Self Placed - Not Authorised - Other Person/s",
        "Self Placed - Not Authorised - Parent/s",
        "Absent - Location Unknown") ~ 1,
      PriorityPlacementExitReason == "CYP Missing" ~ 1,
      PriorityPlacementType == "Absent - Location Unknown" ~ 1,
      ActualPlacementExitReason == "CYP Has Self Restored" ~ 1,
      TRUE ~ 0),
    # create flag for independent living placement
    independent_living_placement = case_when(
      ActualPlacementType %in% c(
        "Independent living") ~ 1,
      TRUE ~ 0
    ),
    # create flag for self placed and/or runaway at any point after age 16
    self_placed_missing_absent_from_placement_after_16 = case_when(
      (self_placed_missing_absent_from_placement == 1 &
         int_overlaps(
           interval(
             date_16 - days(1),
             date_18),
           interval(
             placement_start_date,
             placement_end_date_adjusted
           )
         )) ~ 1,
      TRUE ~ 0)
    ) 

saveRDS(oohc_data_summarise_step_four, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_four.RDS")
#oohc_data_summarise_step_four <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_four.RDS")

# step five: identify distinct care categories
oohc_data_summarise_step_five <- oohc_data_summarise_step_four |>
  # treat every individual separately
  group_by(
    childstory_id
  ) |>
  # create a flag for each placement
  mutate(
    placement_id = row_number()
  ) |>
  ungroup() |>
  group_by(
    childstory_id,
    care_category_start_date,
    care_category_end_date
  ) |>
  summarise(
    placement_id = min(placement_id),
    .groups = "drop"
  ) |>
  mutate(
    # create flag for distinct care category
    distinct_care_category_flag = 1
  ) |>
  ungroup() |>
  select(
    childstory_id,
    placement_id,
    distinct_care_category_flag
  )

saveRDS(oohc_data_summarise_step_five, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_five.RDS")
#oohc_data_summarise_step_five <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_five.RDS")

# step six: merge unique care categories back into master data
oohc_data_summarise_step_six <- oohc_data_summarise_step_four |>
  group_by(
    childstory_id
  ) |>
  # create placement id
  mutate(
    placement_id = row_number()
  ) |>
  # merge care categories data
  left_join(
    oohc_data_summarise_step_five,
    by = c(
      "childstory_id",
      "placement_id"
    )
  ) |>
  mutate(
    # identify if care category closed before individual is eligible
    closed_care_category_start = case_when(
      distinct_care_category_flag == 1 &
      care_category_end_date <= eligible_age_window_start ~ 1,
      TRUE ~ 0
    ),
    # identify if care care category is open when individual is eligible
    open_care_category_start = case_when(
      distinct_care_category_flag == 1 &
      care_category_end_date >= eligible_age_window_start ~ 1,
      TRUE ~ 0
    ),
    # identify if care category closed before eligible window ended
    closed_care_category_end = case_when(
      distinct_care_category_flag == 1 &
        care_category_end_date < eligible_age_window_end ~ 1,
      TRUE ~ 0
    ),
    # identify if care category is still open when eligible window ended
    open_care_category_end = case_when(
      distinct_care_category_flag == 1 &
        care_category_end_date >= eligible_age_window_start ~ 1,
      TRUE ~ 0
    ),
    # calculate time in care for cases closed when eligible window opens
    care_category_length_months_closed_start = case_when(
      distinct_care_category_flag == 1 &
        closed_care_category_start ~ 
          (interval(
            care_category_start_date, 
            care_category_end_date
            ) %/% months(1)),
      TRUE ~ 0),
    # calculate time in care for cases open
    care_category_length_months_open_start = case_when(
      distinct_care_category_flag == 1 &
        open_care_category_start ~ 
        (interval(
          care_category_start_date, 
          eligible_age_window_start
        ) %/% months(1)),
      TRUE ~ 0),
    # calculate time in care for cases closed when eligible window opens
    care_category_length_months_closed_end = case_when(
      distinct_care_category_flag == 1 &
        closed_care_category_end ~ 
        (interval(
          care_category_start_date, 
          care_category_end_date
        ) %/% months(1)),
      TRUE ~ 0),
    # calculate time in care for cases open
    care_category_length_months_open_end = case_when(
      distinct_care_category_flag == 1 &
        open_care_category_end ~ 
        (interval(
          care_category_start_date, 
          eligible_age_window_end
        ) %/% months(1)),
      TRUE ~ 0),
    # sum time in care at start of eligible window
    time_in_care_start_eligible_window =
      max(care_category_length_months_closed_start) +
      max(care_category_length_months_open_start),
    # sum time in care at end of eligible window
    time_in_care_end_eligible_window = 
      max(care_category_length_months_closed_end) +
      max(care_category_length_months_open_end),
    # flag if eligible at start of window
    eligible_time_in_care_start_flag = case_when(
      time_in_care_start_eligible_window > 12 ~ 1,
      TRUE ~ 0
    ),
    # flag if eligible at end of window
    eligible_time_in_care_end_flag = case_when(
      time_in_care_end_eligible_window > 12 ~ 1,
      TRUE ~ 0
    ),
    # flag if eligible at either end of window
    eligible_time_in_care = case_when(
      eligible_time_in_care_start_flag == 1 |
      eligible_time_in_care_end_flag == 1 ~ 1,
      TRUE ~ 0)
  ) |>
  # drop redundant vars
  select(
    -distinct_care_category_flag,
    -closed_care_category_start,              
    -open_care_category_start,
    -closed_care_category_end,
    -open_care_category_end,
    -care_category_length_months_closed_start,
    -care_category_length_months_open_start,
    -care_category_length_months_closed_end,
    -care_category_length_months_open_end,
    time_in_care_start_eligible_window,    
    time_in_care_end_eligible_window,
    -eligible_time_in_care_start_flag,
    -eligible_time_in_care_end_flag
  )

saveRDS(oohc_data_summarise_step_six, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_six.RDS")
#oohc_data_summarise_step_six <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_six.RDS")

# step seven: identify placement individuals with placement instability
oohc_data_summarise_step_seven <- oohc_data_summarise_step_six |>
  # calculate placement length
  mutate(
    placement_length_days = time_length(
      interval(
        placement_start_date, 
        placement_end_date_adjusted
      ) ,
      "days")
  ) |>
  # treat every individual separately
  group_by(
    childstory_id
  ) |>
  # create a flag for each placement
  mutate(
    placement_id = row_number()
  ) |>
  ungroup() |>
  group_by(
    childstory_id,
    placement_start_date,
    placement_end_date_adjusted,
    placement_length_days
  ) |>
  summarise(
    placement_id = min(placement_id),
    .groups = "drop"
  ) |>
  mutate(
    # create flag for distinct placement > 30 days
    distinct_placement_thirty_days_flag = case_when(
      placement_length_days >= 30 ~ 1,
      TRUE ~ 0)
  ) |>
  ungroup() |>
  # drop redundant columns
  select(
    -placement_start_date,
    -placement_end_date_adjusted
  )

saveRDS(oohc_data_summarise_step_seven, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_seven.RDS")
#oohc_data_summarise_step_seven <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_seven.RDS")

# step eight: identify if eligible for placement instability reasons before or during eligibility window
oohc_data_summarise_step_eight <- oohc_data_summarise_step_six |>
  # join data with unique placement flags
  left_join(
    oohc_data_summarise_step_seven,
    by = c(
      "childstory_id",
      "placement_id"
    )
  ) |>
  mutate(
    eligible_placement_length_lifetime_start_eligibility_window = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          date_of_birth,
          eligible_age_window_start)
      ) & distinct_placement_thirty_days_flag == 1 ~ 1,
      TRUE ~ 0
      ),
    eligible_placement_length_lifetime_end_eligibility_window = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          date_of_birth,
          eligible_age_window_end)
      ) & distinct_placement_thirty_days_flag == 1 ~ 1,
      TRUE ~ 0
      ),
    eligible_placement_length_six_months_before_start_eligibility_window = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          add_with_rollback(eligible_age_window_start, - months(6)),
          eligible_age_window_start)
      ) & distinct_placement_thirty_days_flag == 1 ~ 1,
      TRUE ~ 0
      ),
    eligible_placement_length_twelve_months_before_start_eligibility_window = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          add_with_rollback(eligible_age_window_start, - months(12)),
          eligible_age_window_start)
      ) & distinct_placement_thirty_days_flag == 1 ~ 1,
      TRUE ~ 0
      ),
    eligible_placement_length_during_eligibility_window = case_when(
      int_overlaps(
        interval(
          placement_start_date, 
          placement_end_date_adjusted), 
        interval(
          eligible_age_window_start,
          eligible_age_window_end
          )
      ) & distinct_placement_thirty_days_flag == 1 ~ 1,
      TRUE ~ 0
    ),
    # determine eligibility via lifetime placement instability
    eligible_lifetime_placement_instability = case_when(
      sum(eligible_placement_length_lifetime_start_eligibility_window) > 3 |
      sum(eligible_placement_length_lifetime_end_eligibility_window) > 3  ~ 1,
      TRUE ~ 0
    ),
    # determine eligibility via recent placement instability
    eligible_recent_placement_instability = case_when(
      sum(eligible_placement_length_six_months_before_start_eligibility_window) > 0 |
      sum(eligible_placement_length_twelve_months_before_start_eligibility_window) > 0 |
      sum(eligible_placement_length_during_eligibility_window) > 0 ~ 1,
      TRUE ~ 0
    ),
    # determine eligibility via placement instability
    eligible_placement_instability = case_when(
      eligible_lifetime_placement_instability == 1 |
      eligible_recent_placement_instability == 1 
      ~ 1,
      TRUE ~ 0
    )
  ) |>
  # drop redundant columns
  select(
    -distinct_placement_thirty_days_flag,
    -placement_length_days,
    -eligible_placement_length_lifetime_start_eligibility_window,
    -eligible_placement_length_lifetime_end_eligibility_window,
    -eligible_placement_length_six_months_before_start_eligibility_window,
    -eligible_placement_length_twelve_months_before_start_eligibility_window,
    -eligible_placement_length_during_eligibility_window
  )
   
saveRDS(oohc_data_summarise_step_eight, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_eight.RDS")
oohc_data_summarise_step_eight <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_eight.RDS")
 
# step nine: create flag for last placement is self placed and/or runaway at any point
oohc_data_summarise_step_nine <- oohc_data_summarise_step_eight |>
  mutate(
    # count number of placements
    individual_number_of_placements = max(placement_id),
    # identify the last placements
    last_placement_flag = case_when(
      placement_id == individual_number_of_placements ~ 1,
      TRUE ~ 0),
    # identify if child or young person was self placed, missing or absent from their last placement
    last_placement_self_placed_missing_absent = case_when(
      last_placement_flag == 1 &
      self_placed_missing_absent_from_placement == 1 ~ 1,
      TRUE ~ 0),
    # identify if child or young person was in independent living arrangements at their last placement
    last_placement_independent_living = case_when(
      last_placement_flag == 1 &
      independent_living_placement == 1 ~ 1,
      TRUE ~ 0)
    ) |>
  select(
    childstory_id,
    placement_id,
    last_placement_self_placed_missing_absent,
    last_placement_independent_living
  )

saveRDS(oohc_data_summarise_step_nine, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_nine.RDS")
oohc_data_summarise_step_nine <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_nine.RDS")

# step ten: bring last placement self placed, missing or absent back into master data and change end dates to account for placements that have been closed while legal orders remain open
oohc_data_summarise_step_ten <- oohc_data_summarise_step_eight |>
  # merge flag back into data
  left_join(
    oohc_data_summarise_step_nine,
    by = c(
      "childstory_id",
      "placement_id"
    )
  ) |>
  mutate(
    # identify where legal orders still open at the end the last placement
    prm_legal_orders_open_end_last_placement = case_when(
      permanent_responsibility_minister == 1 &
        #(last_placement_self_placed_missing_absent == 1 | 
        #last_placement_independent_living == 1) &
          (legal_orders_end_date > placement_end_date_adjusted |
          legal_orders_end_date > care_category_end_date) ~ 1,
      TRUE ~ 0
    ),
    # for cyp with last placement self placed, missing or absent change placement_end_date_adjusted to legal orders end date
    placement_end_date_adjusted = case_when(
      #(last_placement_self_placed_missing_absent == 1 | 
       #last_placement_independent_living == 1) &
      prm_legal_orders_open_end_last_placement == 1  ~ legal_orders_end_date,
      TRUE ~ placement_end_date_adjusted
    ),
    # for cyp with last placement self placed, missing or absent change care category end date to legal orders end date
    care_category_end_date_adjusted = case_when(
      #(last_placement_self_placed_missing_absent == 1 | 
      # last_placement_independent_living == 1) &
      prm_legal_orders_open_end_last_placement == 1 ~ legal_orders_end_date,
      TRUE ~ care_category_end_date
    ),
    # adjust placement instability criteria so that someone with a last placement as independent living gets flagged as unstable placement
    eligible_placement_instability = case_when(
      last_placement_independent_living == 1 ~ 1,
      TRUE ~ eligible_placement_instability
    )
  ) 
  
saveRDS(oohc_data_summarise_step_ten, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_ten.RDS")
oohc_data_summarise_step_ten <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_ten.RDS")

# step eleven: determine if individual meets eligibility for pyi during evaluation period
oohc_data_summarise_step_eleven <- oohc_data_summarise_step_ten |>
  group_by(
    childstory_id
  ) |>
  mutate(
    # create flag for eligible for pyi
    meets_pyi_eligibility_criteria = case_when(
      eligible_placement_instability == 1 ~ 1,
      eligible_time_in_care == 1 ~ 1,
      eligible_residential_care_placement == 1 ~ 1,
      eligible_permanent_placement == 1 ~ 1,
      TRUE ~ 0),
    meets_pyi_eligibility_criteria = max(
      meets_pyi_eligibility_criteria)
  ) |>
  mutate(
    # create flag if care category overlaps with evaluation period
    care_category_open_during_evaluation = case_when(
      int_overlaps(
        interval(
          care_category_start_date,
          care_category_end_date_adjusted),
        interval(
          pyi_evaluation_start,
          pyi_evaluation_end)
      ) ~ 1,
      TRUE ~ 0
    ),
    # determine date when individual is eligible for pyi
    eligibility_window_start = case_when(
      care_category_open_during_evaluation == 1 &
        age_eligible_during_evaluation == 1 &
        meets_pyi_eligibility_criteria == 1 ~ pmax(pyi_evaluation_start, eligible_age_window_start, care_category_start_date),
      TRUE ~ NA
    ), 
    # determine date when individual is no longer eligible for pyi
    eligibility_window_end = case_when(
      care_category_open_during_evaluation == 1 &
        age_eligible_during_evaluation == 1 &
        meets_pyi_eligibility_criteria == 1 ~ pmin(pyi_evaluation_end, date_18),
      TRUE ~ NA
    ),
    # calculate days eligible
    days_eligible_pyi = case_when(
      care_category_open_during_evaluation == 1 &
        age_eligible_during_evaluation == 1 &
        meets_pyi_eligibility_criteria == 1 ~ as.numeric(
          interval(
            eligibility_window_start,
            eligibility_window_end) %/% days()),
      TRUE ~ 0
    ),
    # flag if ever eligible for pyi during evaluation
    meets_pyi_eligibility_criteria_during_evaluation = case_when(
      days_eligible_pyi > 0 ~ 1,
      TRUE ~ 0
    )) 
  
saveRDS(oohc_data_summarise_step_eleven, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_eleven.RDS")
oohc_data_summarise_step_eleven <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_eleven.RDS")

# add any other variables that might be of interest for balance checks
oohc_data_summarise_step_twelve <- oohc_data_summarise_step_eleven |>
  mutate(
    # flag if placement ever ended for disruptive behaviour  
    placement_ended_disruptive_behaviour = case_when(
      ActualPlacementExitReason %in% c(
        "Disruption Involving CYP") |
      PriorityPlacementExitReason %in% c(
        "Disruption Involving CYP"
      ) ~ 1,
      TRUE ~ 0
    ),
    # flag if placement ended for disruptive behaviour after 16
    placement_ended_disruptive_behaviour_after_16 = case_when(
      (placement_ended_disruptive_behaviour == 1 &
       int_overlaps(
         interval(
           date_16 - days(1),
           date_18),
         interval(
           placement_start_date,
           placement_end_date_adjusted
         )
       )) ~ 1,
    TRUE ~ 0
    ),
    # flag if placement ever broke down  
    placement_breakdown = case_when(
      ActualPlacementExitReason %in% c(
        "Placement Breakdown") |
        PriorityPlacementExitReason %in% c(
          "Placement Breakdown"
        ) ~ 1,
      TRUE ~ 0
    ),
    # flag if placement broke down after 16
    placement_breakdown_after_16 = case_when(
      (placement_breakdown == 1 &
         int_overlaps(
           interval(
             date_16 - days(1),
             date_18),
           interval(
             placement_start_date,
             placement_end_date_adjusted
           )
         )) ~ 1,
      TRUE ~ 0
    ),
    # flag if placement ever ended due to allegation against carer 
    placement_breakdown_allegation_against_carer = case_when(
      ActualPlacementExitReason %in% c(
        "Allegation Against Carer") |
        PriorityPlacementExitReason %in% c(
          "Allegation Against Carer"
        ) ~ 1,
      TRUE ~ 0
    ),
    # flag if placement ended for disruptive behaviour after 16
    placement_breakdown_allegation_against_carer_after_16 = case_when(
      (placement_breakdown_allegation_against_carer == 1 &
         int_overlaps(
           interval(
             date_16 - days(1),
             date_18),
           interval(
             placement_start_date,
             placement_end_date_adjusted
           )
         )) ~ 1,
      TRUE ~ 0
    ),
    # flag if independent living placement after 16
    independent_living_placement_after_16 = case_when(
      (independent_living_placement == 1 &
         int_overlaps(
           interval(
             date_16 - days(1),
             date_18),
           interval(
             placement_start_date,
             placement_end_date_adjusted
           )
         )) ~ 1,
      TRUE ~ 0
    )
  )

saveRDS(oohc_data_summarise_step_twelve, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_twelve.RDS")
oohc_data_summarise_step_twelve <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_twelve.RDS")

# step thirteen: add location and pyi location-based eligibility flags
oohc_data_summarise_step_thirteen <- oohc_data_summarise_step_twelve |>
  mutate(
    # create district flag 
    district_catchment = case_when(
      District == "Central Coast District" ~ "Hunter & Central Coast",         
      District == "Far West District" ~ "Far West",               
      District == "Hunter and Central Coast Region" ~ "Hunter & Central Coast",
      District == "Hunter District" ~ "Hunter & Central Coast",                 
      District == "Hunter New England District" ~ "Hunter & Central Coast",     
      District == "Illawarra Shoalhaven District" ~ "Illawarra Shoalhaven & Southern NSW",   
      District == "Metro Central Region" ~ "Northern Sydney",            
      District == "Metro South West Region" ~ "South Western Sydney",        
      District == "Mid North Coast District" ~ "Mid North Coast & Northern NSW",        
      District == "Murrumbidgee District" ~ "Murrimbidgee",           
      District == "Nepean Blue Mountains District" ~ "Nepean Blue Mountains",  
      District == "New England District" ~ "New England",            
      District == "Northern NSW District" ~ "Mid North Coast & Northern NSW",          
      District == "Northern Sydney District" ~ "Northern Sydney",        
      District == "South Eastern Sydney District" ~ "South Eastern Sydney",   
      District == "South Western Sydney District" ~ "South Western Sydney",   
      District == "Southern NSW District" ~ "Illawarra Shoalhaven & Southern NSW",           
      District == "Southern Region" ~ "Illawarra Shoalhaven & Southern NSW",                
      District == "Statewide Services District" ~ "Statewide Services",     
      District == "Sydney District" ~ "Sydney",                 
      District == "Western NSW District" ~ "Western NSW",            
      District == "Western Sydney District" ~ "Western Sydney",
      TRUE ~ District),
    # identify location at start of eligible window
    district_catchment_eligible_age_start = case_when(
      eligibility_window_start %within%
        interval(
          placement_start_date,
          placement_end_date_adjusted
        ) ~ district_catchment,
      TRUE ~ NA
    ),
    # create flag for pyi eligible based on location
    pyi_eligible_location = case_when(
      district_catchment == "Hunter & Central Coast" ~ 1,
      district_catchment == "Illawarra Shoalhaven & Southern NSW" ~ 1,
      district_catchment == "Mid North Coast & Northern NSW" ~ 1,      
      district_catchment == "Nepean Blue Mountains" ~ 1,               
      district_catchment == "New England" ~ 1,          
      district_catchment == "South Western Sydney" ~ 1,                
      district_catchment == "Statewide Services" & InternalOrganisationBusinessUnit == "Metro ISS" ~ 1,       
      district_catchment == "Western NSW" ~ 1,
      TRUE ~ 0)
    ) |>
  # fill NA values in district_catchment_eligible_age_start
  fill(
    district_catchment_eligible_age_start,
    .direction = "downup"
  ) |>
  # drop three cases with missing location information
  filter(
      !district_catchment %in% c("", "-")
    ) |>
  # remove redundant location columns
  select(
    -District,
    -InternalOrganisationBusinessUnit
  )

saveRDS(oohc_data_summarise_step_thirteen, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_thirteen.RDS")
oohc_data_summarise_step_thirteen <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_thirteen.RDS")

# step fourteen: summarise data
oohc_data_summarise_step_fourteen <- oohc_data_summarise_step_thirteen |>
  select(
    childstory_id,
    psp_slk,
    slk_trimmed,
    district_catchment,                                 
    district_catchment_eligible_age_start,
    pyi_eligible_location,
    pyi_flag,
    matching_pool_flag,
    date_of_birth,
    date_19,
    date_18,
    date_16,
    male,
    female,
    intersex_other,
    aboriginal,
    residential_care_placement,  
    eligible_residential_care_placement,
    eligible_foster_care_placement,
    eligible_kinship_care_placement,                    
    eligible_permanent_placement,
    permanent_responsibility_minister_eligiblity_period,
    self_placed_missing_absent_from_placement,          
    independent_living_placement, 
    independent_living_placement_after_16,
    self_placed_missing_absent_from_placement_after_16,
    time_in_care_start_eligible_window,
    eligible_time_in_care,
    eligible_lifetime_placement_instability,
    eligible_recent_placement_instability,
    eligible_placement_instability,
    eligibility_window_start,                            
    eligibility_window_end, 
    days_eligible_pyi,
    last_placement_self_placed_missing_absent,
    last_placement_independent_living,
    meets_pyi_eligibility_criteria,
    meets_pyi_eligibility_criteria_during_evaluation,
    placement_ended_disruptive_behaviour,
    placement_ended_disruptive_behaviour_after_16,
    placement_breakdown,
    placement_breakdown_after_16,                      
    placement_breakdown_allegation_against_carer,
    placement_breakdown_allegation_against_carer_after_16
  ) |>
  # group by individual
  group_by(
    childstory_id
  ) |>
  # fill some variables so that they can be summarised
  fill(
    eligibility_window_start, .direction = "downup"
  ) |> 
  fill(                         
    eligibility_window_end, .direction = "downup"
  ) |>
  mutate(
    # where there are conflicts between rows take min date
    eligibility_window_start = min(eligibility_window_start),
    # where there are conflicts between rows take max date
    eligibility_window_end = min(eligibility_window_end)
  ) |>
  # summarise values
  mutate(
    # replace NA's in pyi_flag
    pyi_flag = max(pyi_flag),
    pyi_flag = replace_na(
      pyi_flag, 0
      ),
    # take max value of all flags
    across(
      where(
        is.numeric), max)
  ) |>
  # summarise all by taking distinct values
  distinct()

saveRDS(oohc_data_summarise_step_fourteen, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_fourteen.RDS")
oohc_data_summarise_step_fourteen <- readRDS("P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_fourteen.RDS")

# step fifteen: filter ineligible cases
oohc_data_summarise_step_fifteen <- oohc_data_summarise_step_fourteen |>
  mutate(
    in_pyi_location_did_not_receive = case_when(
      pyi_flag == 0 & pyi_eligible_location == 1 ~ 1,
      TRUE ~ 0
    )
  ) |>
  filter(
    # drop individuals who do not meeting the eligibility criteria during the evaluation
    meets_pyi_eligibility_criteria_during_evaluation == 1,
    # drop individuals who are in areas where pyi is provided but did not receive services
    in_pyi_location_did_not_receive == 0,
    # drop individuals where date_19 > extraction date (2021-06-30)
    !date_19 > ymd("2021-06-30")
  ) |>
  rename(
    aihw_slk = psp_slk
  )

saveRDS(oohc_data_summarise_step_fifteen, "P:/pyi/pyi_update/data/processed_data/oohc_data_summarise_step_fifteen.RDS")

  