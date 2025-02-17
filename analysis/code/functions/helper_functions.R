# helper functions for this project

# function to check if packages are installed, install if required and load in groundhog framework
check_install_groundhog <- function(package_list, snapshot_date) {
  # check if groundhog is installed
  if(!require("groundhog", quietly = TRUE)) {
    install.packages("groundhog")
  }
  library("groundhog")
 
  # ensure vector of packages is a character vector
  package_list <- as.character(
    substitute(
      package_list))[-1]
  
  # check what packages are not installed
  missing_packages <- package_list[!package_list %in% installed.packages()[, "Package"]]
                                   
  # install packages if required
  if(length(missing_packages) > 0) {
    message(
      "Initialising weapons, set phasers to stun, install missing packages: ",
      paste(missing_packages, collapse = ", "))
    for (package in missing_packages) {
      groundhog.library(
        package,
        snapshot_date)
    }
  }
    
  # load packages using groundhog   
  for (package in missing_packages) {
      groundhog.library(
        package,
        snapshot_date
      )
    }
  }

# set timezone
Sys.setenv(TZ = "Australia/Sydney")

# set up workspace so that I can call common tidyverse functions
workspace_setup <- function() {
  # import commonly used functions
  import_list <- list(
    # dplyr
    mutate = dplyr::mutate,
    select = dplyr::select,
    summarise = dplyr::summarise,
    rename = dplyr::rename,
    arrange = dplyr::arrange,
    filter = dplyr::filter,
    case_when = dplyr::case_when,
    distinct = dplyr::distinct,
    full_join = dplyr::full_join,
    left_join = dplyr::left_join,
    group_by = dplyr::group_by,
    across = dplyr::across,
    all_of = dplyr::all_of,
    n = dplyr::n,
    slice = dplyr::slice,
    n_distinct = dplyr::n_distinct,
    row_number = dplyr::row_number,
    ungroup = dplyr::ungroup,
    # stringr
    str_sub = stringr::str_sub,
    str_starts = stringr::str_starts,
    # lubridate
    ymd = lubridate::ymd,
    mdy = lubridate::mdy,
    years = lubridate::years,
    month = lubridate::month,
    days = lubridate::days,
    interval = lubridate::interval,
    int_overlaps = lubridate::int_overlaps,
    `%within%` = lubridate::`%within%`,
    time_length = lubridate::time_length,
    add_with_rollback = lubridate::add_with_rollback,
    # purrr
    map = purrr::map,
    map2 = purrr::map2,
    # tidyr
    replace_na = tidyr::replace_na,
    unnest = tidyr::unnest,
    fill = tidyr::fill,
    # ggplot2
    ggplot = ggplot2::ggplot,
    aes = ggplot2::aes,
    geom_histogram = ggplot2::geom_histogram
  )
  
  # attach functions to global environment
  list2env(
    import_list,
    envir = .GlobalEnv
  )
}

# function to check if all values in an id grouping are similar e.g., for sex and aboriginality
check_record_consistency <- function(data, id_column, category_column) {
  
  # find ids with inconsistent categories
  inconsistent_identifiers <- data |>
    group_by(
      across(all_of(id_column))
      ) |>
    summarise(
      count_distinct_categories = n_distinct(across(all_of(category_column))),
      categories = paste(unique(across(all_of(category_column))), collapse = ", "),
      .groups = "drop"
      ) |>
    filter(
      count_distinct_categories > 1
    )
  
  return(inconsistent_identifiers)
}

# function to load series of RDS files
sherlock_reads_files <- function(folder_path, common_prefix) {
  # list all files in directory
  files <- list.files(
    path = folder_path,
    pattern = paste0("^", common_prefix, ".*\\.RDS$"),
    full.names = TRUE
  )
  
  # list sherlock inspired progress messages
  progress_messages <- c(
    "The game is afoot! Reading your files...",
    "No better way to solve a mystery than with proper data",
    "These files hold the key to our investigation",
    "Observe carefully Watson",
    "These data will help us eliminate the impossible"
  )
  
  # create progress bar
  progress <- progress::progress_bar$new(
    format = " :message [:bar] :percent",
    total = length(files),
    width = 60,
    clear = FALSE
  )
  
  # read files with progress updates
  for(file in files) {
    progress$tick(tokens = list(
      message = sample(progress_messages, 1)
    ))
  # name data frame
  object_name <- tools::file_path_sans_ext(basename(file))
  
  # read files into global environment
  assign(object_name, readRDS(file), envir = .GlobalEnv)
  }
}

# function to load series of RDS files
sherlock_reads_files_ending <- function(folder_path, common_suffix) {
  # list all files in directory
  files <- list.files(
    path = folder_path,
    pattern = paste0(".*", common_suffix, ".*\\.RDS$"),
    full.names = TRUE
  )
  
  # list sherlock inspired progress messages
  progress_messages <- c(
    "The game is afoot! Reading your files...",
    "No better way to solve a mystery than with proper data",
    "These files hold the key to our investigation",
    "Observe carefully Watson",
    "These data will help us eliminate the impossible"
  )
  
  # create progress bar
  progress <- progress::progress_bar$new(
    format = " :message [:bar] :percent",
    total = length(files),
    width = 60,
    clear = FALSE
  )
  
  # read files with progress updates
  for(file in files) {
    progress$tick(tokens = list(
      message = sample(progress_messages, 1)
    ))
    # name data frame
    object_name <- tools::file_path_sans_ext(basename(file))
    
    # read files into global environment
    assign(object_name, readRDS(file), envir = .GlobalEnv)
  }
}
  
 