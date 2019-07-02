#' Determine activity according to product formed
#'
#' \code{parse_activity} will determine specific and non-specific activity, expressed in terms of product formed,
#' uisng a capillary migration model calculated by \code{fragr::model_migration}
#'
#' The regression model variables are used to create a predictive function,
#' which allows \code{parse_activity} to select the peak closest to the theoretical or predicted size
#' and label this peak as the specific or expected activity in relation to all activity.
#' The "act" variable represents the activity determined from a peak closest to the theoretical predicted size,
#' which is determined by the regression model variables and alllowed deviation reg_limit arguments.
#' The "act_bp" variable identifies the called size of the peak which was used to generate the "act" variable.
#' The "offsense" variable represents the sum of all other activity besides "act" on the same strand.
#' The "offanti" variable is only included in situations where the substate is double-stranded,
#' and represents the sum of all activity falling on the opposite strand as "act" and "offsense".
#'
#' \code{parse_activity}
#'
#' @seealso [fragr::model_migration]
#'
#' @param nested_df a nested data frame of CE data with all peak data lying in the list column 'data'
#' @param reg_vars a named vector of numeric values from \code{lm} capillary migration model,
#' "slope" and "intercept" must be present
#' @param reg_limit a numeric value defining how far the product peak may deviate from the model
#' @param substrate_cutoff a numeric value indicating the bp size above which all peaks should represent product
#' @param top_n a numeric value for noise reduction selecting the top number of peaks to analyze by peak area,
#' defaults to five (5) unless otherwise specified
#' @param window_size a numeric value indicating the ±window around the specific peak where off-target can be called on the same strand,
#' defaults to two (2) bp above and below "act_bp" unless otherwise specified
#' @return returns the nested data frame with added variables "act", "act_bp", "offsense",
#' and "offanti"
#' @export
parse_activity <- function(nested_df, reg_vars, reg_limit, substrate_cutoff, top_n = 5, window_size = 2) {

  # add activity variables
  if ("ds" %in% nested_df$ssds) {
    nested_df <- nested_df %>%
      tibble::add_column(
        act = as.numeric(NA),
        offsense = as.numeric(NA),
        offanti = as.numeric(NA)
      )
  } else {
    nested_df <- nested_df %>%
      tibble::add_column(
        act = as.numeric(NA),
        offsense = as.numeric(NA)
      )
  }

  # act_bp is being set to NA when there is no specific product peak,
  # I need to have a way to set act_bp to be the model hypothetical act_bp in the event that there is no specific activity,
  # or have a way for offsense to be predicted around that region

  ### SPECIFIC
  # function to extract a specific activity values from the list column "data"
  extract_act <- function(df, dye) {
    dye <- rlang::enquo(dye)
    df %>%
      dplyr::filter(dye == !!dye) %>%
      dplyr::select(rel_area, bp) %>%
      dplyr::top_n(top_n, rel_area) %>%
      dplyr::filter(
        dplyr::between(
          bp,
          ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]]) - reg_limit,
          ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]]) + reg_limit
        ) | bp > substrate_cutoff
      ) %>%
      dplyr::filter(
        abs(bp - ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]])) ==
          min(abs(bp - ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]])))
      ) %>%
      dplyr::pull(var = rel_area) %>%
      (function(v) {
        dplyr::if_else(
          length(v) == 0,
          as.numeric(NA),
          round(mean(v), digits = 2),
          missing = as.numeric(NA)
        )
      })
  }

  # function to extract a specific activity values from the list column "data"
  extract_act_bp <- function(df, dye) {
    dye <- rlang::enquo(dye)
    df %>%
      dplyr::filter(dye == !!dye) %>%
      dplyr::select(rel_area, bp) %>%
      dplyr::top_n(top_n, rel_area) %>%
      dplyr::filter(
        dplyr::between(
          bp,
          ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]]) - reg_limit,
          ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]]) + reg_limit
        ) | bp > substrate_cutoff
      ) %>%
      dplyr::filter(
        abs(bp - ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]])) ==
          min(abs(bp - ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]])))
      ) %>%
      dplyr::pull(var = bp) %>%
      (function(v) {
        dplyr::if_else(
          length(v) == 0,
          ((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]]),
          round(mean(v), digits = 2),
          missing = as.numeric(NA)
        )
      })
  }

  # calculate specific activity
  for (i in seq_along(nested_df[["FWRV"]])) {
    if (identical(nested_df[["FWRV"]][[i]], "FW")) {
      nested_df[["act"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_act("B")
      nested_df[["act_bp"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_act_bp("B")
    } else if (identical(nested_df[["FWRV"]][[i]], "RV") && identical(nested_df[["ssds"]][[i]], "ds")) {
      nested_df[["act"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_act("G")
      nested_df[["act_bp"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_act_bp("G")
    } else if (identical(nested_df[["FWRV"]][[i]], "RV") && identical(nested_df[["ssds"]][[i]], "ss")) {
      nested_df[["act"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_act("B")
      nested_df[["act_bp"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_act_bp("B")
    } else {
      nested_df[["act"]][[i]] <- as.numeric(NA)
      nested_df[["act_bp"]][[i]] <- as.numeric(((reg_vars[["slope"]] * (nested_df[["prod_size"]][[i]])) + reg_vars[["intercept"]]))
    }
  }

  ### OFF-SENSE
  # function to extract a non-specific sense activity values from the list column "data"
  extract_offsense <- function(df, dye) {
    dye <- rlang::enquo(dye)
    df %>%
      dplyr::filter(dye == !!dye) %>%
      dplyr::select(rel_area, bp) %>%
      dplyr::filter(bp < substrate_cutoff) %>%
      dplyr::filter(
        !(dplyr::between(
          bp,
          nested_df[["act_bp"]][[i]] - window_size,
          nested_df[["act_bp"]][[i]] + window_size
        ))
      ) %>%
      dplyr::pull(var = rel_area) %>%
      (function(v) {
        if (length(v) > 0) {
          round(sum(v), digits = 2)
        } else {
          0
        }
      })
  }

  # function to extract a non-specific sense bp values from the list column "data"
  extract_offsense_bp <- function(df, dye) {
    dye <- rlang::enquo(dye)
    df %>%
      dplyr::filter(dye == !!dye) %>%
      dplyr::select(rel_area, bp) %>%
      dplyr::filter(bp < substrate_cutoff) %>%
      dplyr::filter(
        !(dplyr::between(
          bp,
          nested_df[["act_bp"]][[i]] - window_size,
          nested_df[["act_bp"]][[i]] + window_size
        ))
      ) %>%
      dplyr::pull(var = bp) %>%
      (function(v) {
        if (length(v) > 0) {
          paste(v, collapse = ", ")
        } else {
          as.character(NA)
        }
      })
  }

  # calculate off-target activity on the same strand
  for (i in seq_along(nested_df[["FWRV"]])) {
    if (identical(nested_df[["FWRV"]][[i]], "FW")) {
      nested_df[["offsense"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_offsense("B")
      nested_df[["offsense_bp"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_offsense_bp("B")
    } else if (identical(nested_df[["FWRV"]][[i]], "RV") && identical(nested_df[["ssds"]][[i]], "ds")) {
      nested_df[["offsense"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_offsense("G")
      nested_df[["offsense_bp"]][[i]] <- nested_df[["data"]][[i]] %>%
        extract_offsense_bp("G")
    } else {
      nested_df[["offsense"]][[i]] <- as.numeric(NA)
      nested_df[["offsense_bp"]][[i]] <- as.character(NA)
    }
  }

  ### OFF-ANTI
  # function to extract a non-specific anti activity values from the list column "data"
  extract_offanti <- function(df, dye) {
    dye <- rlang::enquo(dye)
    df %>%
      dplyr::filter(dye == !!dye) %>%
      dplyr::select(rel_area, bp) %>%
      dplyr::filter(bp < substrate_cutoff) %>%
      dplyr::pull(var = rel_area) %>%
      (function(v) {
        if (length(v) > 0) {
          round(sum(v), digits = 2)
        } else {
          0
        }
      })
  }

  # function to extract a non-specific anti bp values from the list column "data"
  extract_offanti_bp <- function(df, dye) {
    dye <- rlang::enquo(dye)
    df %>%
      dplyr::filter(dye == !!dye) %>%
      dplyr::select(rel_area, bp) %>%
      dplyr::filter(bp < substrate_cutoff) %>%
      dplyr::pull(var = bp) %>%
      (function(v) {
        if (length(v) > 0) {
          paste(v, collapse = ", ")
        } else {
          as.character(NA)
        }
      })
  }

  if ("ds" %in% nested_df$ssds) {
    for (i in seq_along(nested_df[["ssds"]])) {
      if (identical(nested_df[["ssds"]][[i]], "ds") && identical(nested_df[["FWRV"]][[i]], "FW")) {
        nested_df[["offanti"]][[i]] <- nested_df[["data"]][[i]] %>%
          extract_offanti("G")
        nested_df[["offanti_bp"]][[i]] <- nested_df[["data"]][[i]] %>%
          extract_offanti_bp("G")
      } else if (identical(nested_df[["ssds"]][[i]], "ds") && identical(nested_df[["FWRV"]][[i]], "RV")) {
        nested_df[["offanti"]][[i]] <- nested_df[["data"]][[i]] %>%
          extract_offanti("B")
        nested_df[["offanti_bp"]][[i]] <- nested_df[["data"]][[i]] %>%
          extract_offanti_bp("B")
      }
    }
  }
  return(nested_df)
}
