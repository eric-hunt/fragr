#' Import PeakScanner *combined table* .csv files into R
#'
#' \code{read_PeakScanner}
#'
#' @param directory_path a path to a directory containing the .csv files
#' @param pattern a regex pattern for selecing files in the directory,
#' defaults to reading all .csv files present
#' @param type a character string, either 'legacy' or 'cloud' indicating where
#' the file was generated (desktop app or online app); defaults to 'legacy'
#' @return generates a named list of tibbles,
#' each neamed element is one *combined table* named with its origin file path
#' @export
read_PeakScanner <- function(directory_path, pattern = NULL, type = "legacy") {
  if (!(type %in% c("legacy", "cloud"))) {
    stop("Arg `type` must be either 'legacy' (for desktop app) or 'cloud'.")
  }

  # PeakScanner default columns are not named well for R.
  # These variables provide a way to generate easier to use names.
  legacy_col_info <- list(
    col_types = readr::cols_only(
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
    ),
    col_switch = c(
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
    ),
    col_drop = c(
      "status",
      "type",
      "comment",
      "edit"
    ),
    col_order = c(
      "ABIF_name",
      "sample_id",
      "standard",
      "method",
      "quality",
      "offscale",
      "dye_peak",
      "height",
      "bp",
      "bp_begin",
      "bp_end",
      "bp_width",
      "bp_area",
      "scan",
      "scan_begin",
      "scan_end",
      "scan_width",
      "scan_area"
    )
  )

  col_types <- legacy_col_info$col_types
  col_switch <- legacy_col_info$col_switch
  col_order <- legacy_col_info$col_order

  cloud_col_info <- list(
    col_types = readr::cols_only(
      "Sample File Name" = "c",
      "Dye Color" = "c",
      "Dye , Sample Peak" = "c",
      "Size" = "d",
      "Height" = "d",
      "Area(Data Point)" = "d",
      "Area(Base Pairs)" = "d",
      "Data Point" = "d",
      "Begin Point (Data Point)" = "d",
      "End Point (Data Point)" = "d",
      "Begin Point (Base Pairs)" = "d",
      "End Point (Base Pairs)" = "d"
    ),
    col_switch = c(
      ABIF_name = "Sample File Name",
      dye_color = "Dye Color",
      dye_peak = "Dye , Sample Peak",
      bp = "Size",
      height = "Height",
      scan_area = "Area(Data Point)",
      bp_area = "Area(Base Pairs)",
      scan = "Data Point",
      scan_begin = "Begin Point (Data Point)",
      scan_end = "End Point (Data Point)",
      bp_begin = "Begin Point (Base Pairs)",
      bp_end = "End Point (Base Pairs)"
    ),
    col_drop = c(
      "dye_color"
    ),
    col_order = c(
      "ABIF_name",
      "dye_peak",
      "height",
      "bp",
      "bp_begin",
      "bp_end",
      "bp_area",
      "scan",
      "scan_begin",
      "scan_end",
      "scan_area"
    )
  )

  if (type == "cloud") {
    col_types <- cloud_col_info$col_types
    col_switch <- cloud_col_info$col_switch
    col_order <- cloud_col_info$col_order
  }

  # If no pattern is provided, the funtion will read all .csv files in
  # the directory provided.
  if (missing(pattern) | is.null(pattern)) {
    pattern <- c("\\.csv")
  }

  # Create a list of files from the user-provided directory,
  # by pattern if provided.
  file_list <- list.files(
    path = directory_path,
    pattern = pattern,
    full.names = TRUE
  ) %>%
    purrr::set_names()
  # `purrr` must be specified because `magrittr` also has a set_names
  # which requires two args..
  # This imports each file as an element in a list, then renames the columns,
  # and separates the dye/peak column
  files <- purrr::map(
    file_list,
    readr::read_csv,
    col_types = col_types
  ) %>%
    purrr::map(dplyr::rename, !!!col_switch) %>%
    purrr::map(dplyr::select, col_order) %>%
    purrr::map(
      tidyr::separate, "dye_peak",
      into = c("dye", "peak"),
      sep = "(,)|(, )",
      convert = TRUE
    ) %>%
    purrr::map(tidyr::replace_na, list(bp = 0)) %>%
    purrr::imap(
      ~ tibble::add_column(
        .x,
        csv_name = stringr::str_extract(.y, "[^/]+$"),
        .before = "ABIF_name"
      )
    )

  return(files)
}

#' Import single PeakScanner *combined table* .csv file into R
#'
#' \code{read_PeakScanner_file}
#'
#' @param file_path a path to a directory containing the .csv files
#' @param type a character string, either 'legacy' or 'cloud' indicating where
#' the file was generated (desktop app or online app); defaults to 'legacy'
#' @return generates a named list of tibbles,
#' each neamed element is one *combined table* named with its origin file path
#' @export
read_PeakScanner_file <- function(file_path, type = "legacy") {
  if (!(type %in% c("legacy", "cloud"))) {
    stop("Arg `type` must be either 'legacy' (for desktop app) or 'cloud'.")
  }

  # PeakScanner default columns are not named well for R.
  # These variables provide a way to generate easier to use names.
  legacy_col_info <- list(
    col_types = readr::cols_only(
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
    ),
    col_switch = c(
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
    ),
    col_drop = c(
      "status",
      "type",
      "comment",
      "edit"
    ),
    col_order = c(
      "ABIF_name",
      "sample_id",
      "standard",
      "method",
      "quality",
      "offscale",
      "dye_peak",
      "height",
      "bp",
      "bp_begin",
      "bp_end",
      "bp_width",
      "bp_area",
      "scan",
      "scan_begin",
      "scan_end",
      "scan_width",
      "scan_area"
    )
  )

  col_types <- legacy_col_info$col_types
  col_switch <- legacy_col_info$col_switch
  col_order <- legacy_col_info$col_order

  cloud_col_info <- list(
    col_types = readr::cols_only(
      "Sample File Name" = "c",
      "Dye Color" = "c",
      "Dye , Sample Peak" = "c",
      "Size" = "d",
      "Height" = "d",
      "Area(Data Point)" = "d",
      "Area(Base Pairs)" = "d",
      "Data Point" = "d",
      "Begin Point (Data Point)" = "d",
      "End Point (Data Point)" = "d",
      "Begin Point (Base Pairs)" = "d",
      "End Point (Base Pairs)" = "d"
    ),
    col_switch = c(
      ABIF_name = "Sample File Name",
      dye_color = "Dye Color",
      dye_peak = "Dye , Sample Peak",
      bp = "Size",
      height = "Height",
      scan_area = "Area(Data Point)",
      bp_area = "Area(Base Pairs)",
      scan = "Data Point",
      scan_begin = "Begin Point (Data Point)",
      scan_end = "End Point (Data Point)",
      bp_begin = "Begin Point (Base Pairs)",
      bp_end = "End Point (Base Pairs)"
    ),
    col_drop = c(
      "dye_color"
    ),
    col_order = c(
      "ABIF_name",
      "dye_peak",
      "height",
      "bp",
      "bp_begin",
      "bp_end",
      "bp_area",
      "scan",
      "scan_begin",
      "scan_end",
      "scan_area"
    )
  )

  if (type == "cloud") {
    col_types <- cloud_col_info$col_types
    col_switch <- cloud_col_info$col_switch
    col_order <- cloud_col_info$col_order
  }

  # This imports one single file, then renames the columns,
  # and separates the dye/peak column
  file <- readr::read_csv(file_path, col_types = col_types) %>%
    dplyr::rename(!!!col_switch) %>%
    dplyr::select(col_order) %>%
    tidyr::separate(
      "dye_peak",
      into = c("dye", "peak"),
      sep = "(,)|(, )",
      convert = TRUE
    ) %>%
    tidyr::replace_na(list(bp = 0)) %>%
    tibble::add_column(
      csv_name = stringr::str_extract(file_path, "[^/]+$"),
      .before = "ABIF_name"
    )

  return(file)
}
