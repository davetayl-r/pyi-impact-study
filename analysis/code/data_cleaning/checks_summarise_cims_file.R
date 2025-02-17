# make sense of cims file 

# load helper functions
source("./code/functions/helper_functions.R")

# load required packages using custom helper function
required_packages <- c(
  "tidyverse",
  "readxl",
  "stringr"
)

snapshot_date <- "2024-11-01"

check_install_groundhog(
  required_packages,
  snapshot_date
)

# set up workspace preferences using custom helper function
workspace_setup()

# check frequencies for different counting rules
cims_data_summarise_step_two |>
  select(
    requires_housing_assistance_original,
    requires_housing_assistance,
    unsheltered_homelessness_current,
    unsheltered_homelessness_last_week,
    Previously_Homeless_Mth1,
    Previously_Homeless_Mth2,
    Short_Term_Accom_ind
  ) |>
  ungroup() |>
  summarise(
    requires_housing_assistance_original = sum(requires_housing_assistance_original),
    requires_housing_assistance = sum(requires_housing_assistance),
    unsheltered_homelessness_current = sum(unsheltered_homelessness_current),
    unsheltered_homelessness_last_week = sum(unsheltered_homelessness_last_week),
    Previously_Homeless_Mth1 = sum(Previously_Homeless_Mth1),
    Previously_Homeless_Mth2 = sum(Previously_Homeless_Mth2),
    Short_Term_Accom_ind = sum(Short_Term_Accom_ind)
  )
