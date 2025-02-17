# the purpose of this code is to check if we can locate the kids who received PYI in the oohc data from the PSP evaluation

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

mapped_ids <- mapped_ids_data |>
  # change names for consistency
  rename(
    pyi_slk = SLK,
    childstory_id = CS.ID
  ) |>
  # drop redundant vars
  select(
    pyi_slk,
    childstory_id
  ) |>
  # remove final letter of slk
  mutate(
    slk_trimmed = case_when(
      !is.na(pyi_slk) ~ str_sub(pyi_slk, 1, 13),
      TRUE ~ NA_character_
    )
  ) |>
  # deduplicate
  distinct()

# load pyi matching pool and extract ids
pyi_count_data_location <- "P:/pyi/pyi_update/data/raw_data/pyi_evaluation_matching_pool.rds"
pyi_count_data <- readRDS(pyi_count_data_location)

matching_pool <- pyi_count_data |>
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
    matching_pool = 1
  ) |>
  # add childstory id from eligible pool
  mutate(
    childstory_id = case_when(
      slk_trimmed %in% mapped_ids$slk_trimmed ~ mapped_ids$childstory_id[match(slk_trimmed, mapped_ids$slk_trimmed)],
      TRUE ~ NA_character_)
  )

# load oohc oohc data
oohc_data_location <- "P:/PYI/pyi_update/data/raw_data/psp_oohc_final.csv"
oohc_data <- read.csv(oohc_data_location) 

oohc_ids <- oohc_data |>
  rename(
    oohc_slk = SLK,
    childstory_id = ChildStoryID
  ) |>
  select(
    oohc_slk,
    childstory_id
  ) |>
  distinct(
    oohc_slk,
    childstory_id
  ) |>
  # remove final letter of slk
  mutate(
    slk_trimmed = case_when(
      !is.na(oohc_slk) ~ str_sub(oohc_slk, 1, 13),
      TRUE ~ NA_character_
    )
  )

# identify if there are any missing pyi records in the oohc oohc file
# missing n=3 records
missing_records <- matching_pool |>
  filter(
    pyi_flag == 1,
    !slk_trimmed %in% oohc_ids$slk_trimmed &
      !childstory_id %in% oohc_ids$childstory_id
  ) |>
  select(
    pyi_flag,
    slk_trimmed,
    childstory_id
  )
missing_records
