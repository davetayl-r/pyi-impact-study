# this code joins the pyi matching pool to the psp data

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

# read in list of linked SLKs and childstory ids
mapped_ids_data_location <- "P:/pyi/pyi_update/data/raw_data/slk_childstory_id_mapped.csv"
mapped_ids_data <- read.csv(mapped_ids_data_location)

# read matching pool from pyi project
matching_pool_data_location <- "P:/pyi/PYI/analysis_data/flat_files/Propensity_Score_Match_OOHC_final.rds"
matching_pool_data <- readRDS(matching_pool_data_location)

# read oohc file 
oohc_data_location <- "P:/psp/PSP/Final PSP Raw Data/OOHC Final.csv"
oohc_data <- read.csv(oohc_data_location) 

# extract pyi_slk and childstory ids from those who were ever eligible for pyi
mapped_ids <- mapped_ids_data |>
  rename(
    pyi_slk = SLK,
    childstory_id = CS.ID
  ) |>
  select(
    pyi_slk,
    childstory_id
  ) |>
  # force NA for childstory ids that do not confirm to desired structure
  mutate(
    childstory_id = case_when(
      !str_starts(childstory_id, "C-") ~ NA_character_,
      TRUE ~ childstory_id
    )
  ) |>
  # remove final letter of slk
  mutate(
    slk_trimmed = case_when(
      !is.na(pyi_slk) ~ str_sub(pyi_slk, 1, 13),
      TRUE ~ NA_character_
    ),
    ever_eligible_pyi = 1
  ) |>
  # drop pyi slk
  select(
    -pyi_slk
  )

# load matching pool and extract slk and pyi flag
linked_matching_pool <- matching_pool_data |>
  select(
    Statistical.Linkage.key,
    PYI_kid_flag
  ) |>
  rename(
    pyi_slk = Statistical.Linkage.key,
    pyi_flag = PYI_kid_flag
  ) |>
  mutate(
    slk_trimmed = case_when(
      !is.na(pyi_slk) ~ str_sub(pyi_slk, 1, 13),
      TRUE ~ NA_character_
    ),
    matching_pool_flag = 1
  ) |>
  # drop pyi slk
  select(
    -pyi_slk
  ) |>
  # add childstory id from eligible pool
  mutate(
    childstory_id = case_when(
      slk_trimmed %in% mapped_ids$slk_trimmed ~ mapped_ids$childstory_id[match(slk_trimmed, mapped_ids$slk_trimmed)],
      TRUE ~ NA_character_)
  )

# check that childstory id's were mapped correctly - we expect 269
linked_matching_pool |>
  summarise(
    check_childstory_ids_equals_269 = sum(!is.na(childstory_id))
  )

# check pyi flags equal 297
linked_matching_pool |>
  summarise(
    check_pyi_count_equals_297 = sum(pyi_flag, na.rm = TRUE)
  )

# link pyi data to oohc file
linked_oohc_file <- oohc_data |>
  rename(
    psp_slk = SLK,
    childstory_id = ChildStoryID
  ) |>
  # remove final letter of slk
  mutate(
    slk_trimmed = case_when(
      !is.na(psp_slk) ~ str_sub(psp_slk, 1, 13),
      TRUE ~ NA_character_)
  ) |>
  # import pyi flag
  mutate(
    pyi_flag = case_when(
      slk_trimmed %in% linked_matching_pool$slk_trimmed[linked_matching_pool$pyi_flag == 1] ~ 1,
      childstory_id %in% linked_matching_pool$childstory_id[linked_matching_pool$pyi_flag == 1] ~ 1,
      slk_trimmed %in% linked_matching_pool$slk_trimmed[linked_matching_pool$pyi_flag == 0] ~ 0,
      childstory_id %in% linked_matching_pool$childstory_id[linked_matching_pool$pyi_flag == 0] ~ 0
    )
  ) |>
  # import matching pool flag
  mutate(
    matching_pool_flag = case_when(
      slk_trimmed %in% linked_matching_pool$slk_trimmed[linked_matching_pool$matching_pool_flag == 1] ~ 1,
      childstory_id %in% linked_matching_pool$childstory_id[linked_matching_pool$matching_pool_flag == 1] ~ 1
    )
  ) |>
  # drop redundant vars
  select(
    -financial_created_date,
    -payment_amount_ex_gst,
    -payment_gst_amount,
    -payment_amount_incl_gst,
    -service_type_2, 
    -service_type_3_service,
    -date_paid)

# check to see if all 295 pyi kids are in this file - dim() should equal [1] 295 3
linked_oohc_file |>
  select(
    slk_trimmed,
    childstory_id,
    pyi_flag
  ) |>
  distinct() |>
  filter(pyi_flag == 1) |>
  dim()

# export linked oohc file
saveRDS(
  linked_oohc_file,
  "P:/pyi/pyi_update/data/processed_data/oohc_file_pyi_flags.RDS"
)
