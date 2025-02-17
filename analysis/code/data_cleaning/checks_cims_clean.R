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

# linked oohc and pyi file
oohc_file_location <- "P:/pyi/pyi_update/data/processed_data/oohc_file_pyi_flags.RDS"
oohc_file <- readRDS(oohc_file_location)

# cims file location pyi data
pyi_cims_file_location_1 <- "P:/pyi/PYI/analysis_data/flat_files/SHS_CIMS_all_MASTER.rds"
pyi_cims_file_location_2 <- "P:/pyi/PYI/analysis_data/flat_files/SHS_CIMS_all_MASTER_CLEAN.rds"

pyi_cims_1_raw <- readRDS(pyi_cims_file_location_1)
pyi_cims_2_raw <- readRDS(pyi_cims_file_location_2)

# cims file location psp data
psp_cims_file_location_1 <- "P:/psp/PSP/flat_files_final/CIMS/SHS_CIMS_all_MASTER.rds"
psp_cims_file_location_2 <- "P:/psp/PSP/flat_files_final/CIMS/SHS_CIMS_all_MASTER_CLEAN.rds"

psp_cims_1_raw <- readRDS(psp_cims_file_location_1)
psp_cims_2_raw <- readRDS(psp_cims_file_location_2)

# find min and max date ranges for PYI CIMS data 

dim(pyi_cims_1_raw)
dim(pyi_cims_2_raw)

min(pyi_cims_2_raw$spell_start_date) # "2015-03-26"
max(pyi_cims_2_raw$spell_start_date) # "2020-03-31"

# find min and max date range for PSP CIMS data

dim(psp_cims_1_raw)
dim(psp_cims_2_raw)

min(psp_cims_2_raw$spell_start_date) # "2015-01-12"
max(psp_cims_2_raw$spell_start_date) # "2021-06-30"

# get unique slks for pyi from oohc file
pyi_slks <- oohc_file |>
  select(
    pyi_flag,
    matching_pool_flag,
    slk_trimmed
  ) |>
  filter(
    #matching_pool_flag == 1
    pyi_flag == 1
  ) |>
  distinct(
  )

sum(pyi_slks$pyi_flag, na.rm = TRUE)
sum(pyi_slks$matching_pool_flag, na.rm = TRUE)

# trim slks
pyi_cims_1_raw_cleam <- pyi_cims_1_raw |>
#pyi_cims_2_raw_cleam <- pyi_cims_2_raw |>
  mutate(
  slk_trimmed = case_when(
    !is.na(aihw_slk) ~ str_sub(aihw_slk, 1, 13),
    TRUE ~ NA_character_
  )
) 

psp_cims_1_raw_clean <- psp_cims_1_raw |>
#psp_cims_2_raw_clean <- psp_cims_2_raw |>
  mutate(
    slk_trimmed = case_when(
      !is.na(aihw_slk) ~ str_sub(aihw_slk, 1, 13),
      TRUE ~ NA_character_
    )
  ) 

# find number of pyi flag in pyi cims
pyi_slk_original_extract <- pyi_cims_1_raw_cleam |>
#pyi_slk_original_extract <- pyi_cims_2_raw_cleam |>
  select(slk_trimmed) |>
  mutate(slk_check = case_when(
    slk_trimmed %in% pyi_slks$slk_trimmed ~ 1,
    TRUE ~ 0)
  ) |>
  distinct() |>
  summarise(
    pyi_slk_original_extract = sum(slk_check)
  )
pyi_slk_original_extract 

#pyi_slk_original_extract
#1                      113

# get list of slks from pyi oohc file that appear in first cims extract
vector_pyi_slk_original_extract <- pyi_cims_1_raw_cleam |>
#vector_pyi_slk_original_extract <- pyi_cims_2_raw_cleam |>
  select(
    slk_trimmed
  ) |>
  distinct() |>
  mutate(slk_check = case_when(
    slk_trimmed %in% pyi_slks$slk_trimmed ~ 1,
    TRUE ~ 0)
  ) |>
  filter(
    slk_check == 1
  )
# dim 113  

# find overlap of slks from first extract in second extract
extract_overlap_check <- psp_cims_1_raw_clean |>
#extract_overlap_check <- psp_cims_2_raw_clean |>
  select(slk_trimmed) |>
  distinct() |>
  mutate(check = case_when(
    slk_trimmed %in% vector_pyi_slk_original_extract$slk_trimmed ~ 1,
    TRUE ~ 0)
  ) |>
  summarise(
    overlap_check = sum(check)
  )
# 113

# to do
# find 'correct' file
# min/max check date ranges