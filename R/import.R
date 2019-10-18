#' Import PeakScanner *combined table* .csv files into R
#'
#' \code{read_PeakScanner}
#'
#' @param directory_path a path to a directory containing the .csv files
#' @param pattern a regex pattern for selecing files in the directory,
#' defaults to reading all .csv files present
#' @param clean an opinionated removal of some default PeakScanner variables,
#' defaults to TRUE
#' @return generates a named list of tibbles,
#' each neamed element is one *combined table* named with its origin file path
#' @export
read_PeakScanner <- function(directory_path, pattern = NULL, clean = TRUE) {
  if (!(clean %in% c(TRUE, FALSE))) {
    stop("Arg `clean` must be TRUE or FALSE.")
  }

  # PeakScanner default columns are not named well for R.
  # These variables provide a way to generate easier to use names.
  peakScanner_types <- readr::cols_only(
    "Status" = "c",
    "Sample Name" = "c",
    "Sample Type" = "c",
    "Size Standard" = "c",
    "Analysis Method" = "c",
    "Offscale" = "c",
    "Quality" = "c",
    "Dye/Sample Peak" = "c",
    "Sample File Name" = "c",
    "Size" = "d",
    "Height" = "d",
    "Area in Point" = "d",
    "Area in BP" = "d",
    "Data Point" = "d",
    "Begin Point" = "d",
    "Begin BP" = "d",
    "End Point" = "d",
    "End BP" = "d",
    "Width in Point" = "d",
    "Width in BP" = "d",
    "User Comments" = "c",
    "User Edit" = "c"
  )

  peakScanner_colnames <- c(
    Status = "status",
    `Sample Name` = "sample_id",
    `Sample Type` = "type",
    `Size Standard` = "standard",
    `Analysis Method` = "method",
    Offscale = "offscale",
    Quality = "quality",
    `Dye/Sample Peak` = "dye_peak",
    `Sample File Name` = "ABIF_name",
    Size = "bp",
    Height = "height",
    `Area in Point` = "scan_area",
    `Area in BP` = "bp_area",
    `Data Point` = "scan",
    `Begin Point` = "scan_begin",
    `Begin BP` = "bp_begin",
    `End Point` = "scan_end",
    `End BP` = "bp_end",
    `Width in Point` = "scan_width",
    `Width in BP` = "bp_width",
    `User Comments` = "comment",
    `User Edit` = "edit"
  )
  rlang::enquos(peakScanner_colnames)

  peakScanner_rename <- c(
    status = "Status",
    sample_id = "Sample Name",
    type = "Sample Type",
    standard = "Size Standard",
    method = "Analysis Method",
    offscale = "Offscale",
    quality = "Quality",
    dye_peak = "Dye/Sample Peak",
    ABIF_name = "Sample File Name",
    bp = "Size",
    height = "Height",
    scan_area = "Area in Point",
    bp_area = "Area in BP",
    scan = "Data Point",
    scan_begin = "Begin Point",
    bp_begin = "Begin BP",
    scan_end = "End Point",
    bp_end = "End BP",
    scan_width = "Width in Point",
    bp_width = "Width in BP",
    comment = "User Comments",
    edit = "User Edit"
  )
  rlang::enquos(peakScanner_rename)

  # I want to eventually do this step earlier and call the actual names imported from PeakScanner..
  # ..in this way, a user can supply a vector of column names as an argument instead of this opinionated cleaning.
  vars_selected <- c(
    "csv_name", "ABIF_name", "sample_id", "dye", "peak", "height", "bp", "bp_area", "bp_width",
    "scan", "scan_area", "scan_width", "offscale", "quality", "standard", "method"
  )
  rlang::enquos(vars_selected)

  # If no pattern is provided, the funtion will read all .csv files in the directory provided.
  if (missing(pattern) | is.null(pattern)) {
    pattern <- c("\\.csv")
  }

  # Create a list of files from the user-provided directory, by pattern if provided.
  file_list <- list.files(path = directory_path, pattern = pattern, full.names = TRUE) %>%
    purrr::set_names() # `purrr` must be specified because `magrittr` also has a set_names which requires two args..
  # This imports each file as an element in a list, then renames the columns, and separates the dye/peak column
  files <- purrr::map(file_list, readr::read_csv, col_types = peakScanner_types) %>%
    purrr::map(dplyr::rename, !!! peakScanner_rename) %>%
    purrr::map(tidyr::separate, "dye_peak", into = c("dye", "peak"), sep = ", ", convert = TRUE) %>%
    purrr::map(tidyr::replace_na, list(bp = 0)) %>%
    purrr::map2(., names(.), ~ tibble::add_column(.x, csv_name = stringr::str_extract(.y, "(?<=//).*(?=\\.csv)"), .before = "ABIF_name"))

  if (clean) {
    files %>%
      purrr::map(dplyr::arrange, dye, sample_id, peak) %>%
      purrr::map(dplyr::select, c(!!! vars_selected)) %>%
    return()
  } else {
    return(files)
  }
}

#' Import single PeakScanner *combined table* .csv file into R
#'
#' \code{read_PeakScanner_file}
#'
#' @param file_path a path to a directory containing the .csv files
#' @param clean an opinionated removal of some default PeakScanner variables,
#' defaults to TRUE
#' @return generates a named list of tibbles,
#' each neamed element is one *combined table* named with its origin file path
#' @export
read_PeakScanner_file <- function(file_path, clean = TRUE) {
  if (!(clean %in% c(TRUE, FALSE))) {
    stop("Arg `clean` must be TRUE or FALSE.")
  }

  # PeakScanner default columns are not named well for R.
  # These variables provide a way to generate easier to use names.
  peakScanner_types <- readr::cols_only(
    "Status" = "c",
    "Sample Name" = "c",
    "Sample Type" = "c",
    "Size Standard" = "c",
    "Analysis Method" = "c",
    "Offscale" = "c",
    "Quality" = "c",
    "Dye/Sample Peak" = "c",
    "Sample File Name" = "c",
    "Size" = "d",
    "Height" = "d",
    "Area in Point" = "d",
    "Area in BP" = "d",
    "Data Point" = "d",
    "Begin Point" = "d",
    "Begin BP" = "d",
    "End Point" = "d",
    "End BP" = "d",
    "Width in Point" = "d",
    "Width in BP" = "d",
    "User Comments" = "c",
    "User Edit" = "c"
  )

  peakScanner_colnames <- c(
    Status = "status",
    `Sample Name` = "sample_id",
    `Sample Type` = "type",
    `Size Standard` = "standard",
    `Analysis Method` = "method",
    Offscale = "offscale",
    Quality = "quality",
    `Dye/Sample Peak` = "dye_peak",
    `Sample File Name` = "ABIF_name",
    Size = "bp",
    Height = "height",
    `Area in Point` = "scan_area",
    `Area in BP` = "bp_area",
    `Data Point` = "scan",
    `Begin Point` = "scan_begin",
    `Begin BP` = "bp_begin",
    `End Point` = "scan_end",
    `End BP` = "bp_end",
    `Width in Point` = "scan_width",
    `Width in BP` = "bp_width",
    `User Comments` = "comment",
    `User Edit` = "edit"
  )
  rlang::enquos(peakScanner_colnames)

  peakScanner_rename <- c(
    status = "Status",
    sample_id = "Sample Name",
    type = "Sample Type",
    standard = "Size Standard",
    method = "Analysis Method",
    offscale = "Offscale",
    quality = "Quality",
    dye_peak = "Dye/Sample Peak",
    ABIF_name = "Sample File Name",
    bp = "Size",
    height = "Height",
    scan_area = "Area in Point",
    bp_area = "Area in BP",
    scan = "Data Point",
    scan_begin = "Begin Point",
    bp_begin = "Begin BP",
    scan_end = "End Point",
    bp_end = "End BP",
    scan_width = "Width in Point",
    bp_width = "Width in BP",
    comment = "User Comments",
    edit = "User Edit"
  )
  rlang::enquos(peakScanner_rename)

  # I want to eventually do this step earlier and call the actual names imported from PeakScanner..
  # ..in this way, a user can supply a vector of column names as an argument instead of this opinionated cleaning.
  vars_selected <- c(
    "csv_name", "ABIF_name", "sample_id", "dye", "peak", "height", "bp", "bp_area", "bp_width",
    "scan", "scan_area", "scan_width", "offscale", "quality", "standard", "method"
  )
  rlang::enquos(vars_selected)

  # This imports one single file, then renames the columns, and separates the dye/peak column
  file <- readr::read_csv(file_path, col_types = peakScanner_types) %>%
    dplyr::rename(!!! peakScanner_rename) %>%
    tidyr::separate("dye_peak", into = c("dye", "peak"), sep = ", ", convert = TRUE) %>%
    tidyr::replace_na(list(bp = 0)) %>%
    add_column(csv_name = file_path, .before = "ABIF_name")

  if (clean) {
    file %>%
      dplyr::arrange(dye, sample_id, peak) %>%
      dplyr::select(c(!!! vars_selected)) %>%
      return()
  } else {
    return(file)
  }
}
