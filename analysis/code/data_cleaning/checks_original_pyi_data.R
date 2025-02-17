# check original pyi data for missing cases

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

# specify column types
column_types = c(
  "Statistical Linkage key" = "guess",               
  "CYP date of birth" = "date",                     
  "CYP Gender" = "guess",                           
  "Aboriginal Status (FaCS Grouped)" = "guess",      
  "Care Category Start Date" = "date",              
  "Care Category End Date" = "date",               
  "ExitCareType" = "guess",                          
  "District" = "guess",                              
  "Internal Organisation (Business Unit)" = "guess",
  "Priority Placement Purpose" = "guess",            
  "Priority Placement Start Date" = "date",         
  "Priority Placement End Date" = "date",          
  "Priority Placement Exit Reason" = "guess",        
  "Priority Placement Provider - Grouped" = "guess", 
  "Priority Placement Type" = "guess",              
  "Aboriginal status of Carer Household" = "guess",  
  "Parental Responsibility (Ungrouped)" = "guess",   
  "Parental Responsibility (Grouped)" = "guess",    
  "Legal Order Start Date" = "date",                
  "Legal Order End Date" = "date"
)

# read original pyi data
original_pyi_data_location <- "P:/pyi/PYI/raw_data/Childstory_Aug2020/CYP in OOHC from July 2015 onwards 20200819.xlsx"
original_pyi_data <- readxl::read_excel(
  original_pyi_data_location,
  sheet = 2,
  col_types = column_types
  )

# load missing ids from step_five (assumes that this data is loaded)
missing_ids <- step_five |>
  filter(
    pyi_flag == 1,
    !childstory_id %in% ids_in_data$childstory_id
  ) |>
  select(
    slk_trimmed 
  ) |>
  distinct()

# inspect cases individually 
original_pyi_data |>
  filter(
    str_starts(`Statistical Linkage key`, "ICEAK27102001")
  ) |>
  View()

# 1 AN2AS28122000 Last placement recorded as 'system missing'
# 2 ILOAX08092000 Last placement recorded as 'system missing'
# 3 ARUES27112000 Last placement recorded as 'system missing'
# 4 OHSKY16032001 Last placement recorded as 'system missing'
# 5 ARING01072001 Last placement recorded as 'system missing'
# 6 ACYYR11052001 Last placement recorded as 'system missing'
# 7 ONERA22102001 Last placement recorded as 'system missing'
# 8 ICEAK27102001 Young person recorded as incarcerated



names(original_pyi_data)