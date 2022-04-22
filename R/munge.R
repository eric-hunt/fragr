#' Fix instances where substrate is not defined by a single peak
#'
#' \code{munge_subfix}
#'
#' @param nested_df a nested data frame of CE data with all peak data lying in
#' the list column 'data'
#' @return returns the data frame with the product peak(s) merged into one peak
#' @export
munge_subfix <- function(nested_df) {
  nested_df %>%
    dplyr::mutate(
      data = data %>%
        purrr::modify(
          ~ .x %>%
            tidyr::nest(subdata = -c(dye, label)) %>%
            dplyr::mutate(subdata = dplyr::if_else(
              label == "substrate",
              subdata %>% purrr::modify(
                dplyr::mutate,
                rel_area = sum(rel_area)
              ),
              subdata
            )) %>%
            dplyr::mutate(subdata = dplyr::if_else(
              label == "substrate",
              subdata %>% purrr::modify(dplyr::slice_max, height, n = 1),
              subdata
            )) %>%
            dplyr::mutate(subdata = dplyr::if_else(
              label == "substrate",
              subdata %>% purrr::modify(~ .x %>%
                dplyr::mutate(rel_area = dplyr::if_else(
                  rel_area > 1,
                  1.00,
                  rel_area
                ))),
              subdata
            )) %>%
            tidyr::unnest(subdata)
        )
    )
}



#' Put relative measurements in terms of product formed
#'
#' \code{munge_prodify}
#'
#' @param nested_df a nested data frame of CE data with all peak data lying in
#' the list column 'data'
#' @return returns the data frame relative data expressed in terms of
#' percent product formed (i.e. 100% substrate will now be 0% product)
#' @export
# munge_prodify <- function(nested_df, reg_vars, reg_limit, substrate_cutoff) {
munge_prodify <- function(nested_df) {
  for (i in seq_along(nested_df[[1]])) {
    nested_df[["data"]][[i]] <- nested_df[["data"]][[i]] %>%
      dplyr::group_by(dye) %>%
      dplyr::mutate(
        rel_area = dplyr::if_else(
          dplyr::n() == 1 & label == "substrate",
          0.00,
          rel_area
        ),
        label = dplyr::if_else(
          dplyr::n() == 1 & label == "substrate" & rel_area == 0.00,
          "product",
          label
        )
      ) %>%
      dplyr::filter(label == "product") %>%
      dplyr::ungroup()
  }
  return(nested_df)
}
