#' Generate regression model for CE migration
#'
#'
#'
#' \code{parse_migmodel}
#'
#' @param df_list a list of data frames containing CE data imported with \code{fragr::read_PeakScanner}
#' @param channels the one letter channel/dye code (as a character vector) for which to build a migration regression model
#' @param search_var a metadata-containing variable which contains the theoretical product size
#' @param search_pattern a regex pattern for extracting the theoretical product size 'search_var'
#' @param substrate_cutoff a numeric value indicating the bp size above which all peaks should represent product
#' @param limit a numeric value defining how far, in basepairs, the called value an deviate from the theoretical value,
#' defaults to ten (10) unless otherwise specified
#' @return returns a model plot and list of variables..
#' @export
parse_migmodel <- function(df_list, channels, search_var, search_pattern, substrate_cutoff, limit = 10) {
  dye_selection <- rlang::enquo(channels)
  search_var <- rlang::enquo(search_var)
  search_pattern <- rlang::enquo(search_pattern)

  df <- df_list %>%
    dplyr::bind_rows(.id = "file_name") %>%
    tidyr::nest(-c(file_name, csv_name, ABIF_name, sample_id, dye)) %>%
    dplyr::filter(dye %in% !! dye_selection) %>%
    dplyr::mutate(
      bp_actual = as.integer(stringr::str_extract(.data[[rlang::as_name(search_var)]], !! search_pattern))
      # bp_actual = as.integer(stringr::str_extract(.$UQ(search_var), "\\d{2}")) # alternate unquote method
    ) %>%
    dplyr::mutate(
      bp_called = purrr::map_dbl(data, ~ .x %>%
                                   dplyr::top_n(1, height) %>%
                                   dplyr::pull(bp))
    ) %>%
    dplyr::filter(bp_called < substrate_cutoff) %>%
    dplyr::filter(abs(bp_actual - bp_called) < limit)

  model <- df %>%
    lm(bp_called ~ bp_actual, data = .)
  glance <- broom::glance(model)
  tidy <- broom::tidy(model) %>% tibble::column_to_rownames(var = "term")

  message(glue::glue(
    "The following is a summary of the plotted model.
    y is bp_called, x is bp_actual.
    Formula is based on y = mx + b"
  ))
  print(glance)
  print(tidy)

  message(glue::glue(
    "The following list of values can be provided as an argument for other 'fragr' functions.
    It has been assigned to the global environment as 'reg_vars' for convenience."
  ))
  reg_vars <<- list(
    r.squared = round(glance[["r.squared"]], digits = 3),
    intercept = round(tidy["(Intercept)", "estimate"], digits = 3),
    slope = round(tidy["bp_actual", "estimate"], digits = 3)
  ) %>% dput()

  df %>%
    ggplot2::ggplot(ggplot2::aes(x = bp_actual, y = bp_called)) +
    ggplot2::geom_smooth(method = lm, formula = y ~ x) +
    ggplot2::geom_point() +
    ggplot2::annotate(
      geom = "text", x = -Inf, y = Inf, hjust = -1, vjust = 1,
      label = paste(paste0("dye = ", rlang::eval_tidy(dye_selection)),
        "y = mx + b",
        paste0("(R)^2 = ", round(glance[["r.squared"]], digits = 3)),
        sep = "\n"
      )
    ) +
    ggplot2::scale_x_continuous(
      labels = (function(x) round(exp(x))),
      breaks = scales::pretty_breaks()
    )
}
# For testing:
# read_PeakScanner("~/Desktop/Data Consolidation/data_raw/Kevin_comp/noETSSB/data/") %>%
#   model_migration(B, sample_id, "\\d{2}$", 89, 15)
#
# For debugging:
# search_var <- sym(ABIF_name)
# rlang::qq_show(
#   dplyr::mutate(
#     bp_actual = as.integer(stringr::str_extract(.$UQ(search_var), "\\d{2}"))
#   )
# )
