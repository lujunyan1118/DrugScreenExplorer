# Documentation of datasets

#' Example drug screen data on CLL samples
#'
#' This is an example drug screen dataset of CLL samples processed by \code{readScreen()} function.
#' Normalization was not performed.
#'
#' @docType data
#' @usage data(screenData)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(screenData)
"screenData"


#' Example drug screen data on CLL samples with normalization
#'
#' This is an example drug screen dataset of CLL samples processed by \code{readScreen()} function.
#' Normalization was performed based on negative controls.
#'
#' @docType data
#' @usage data(screenData_normalized)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(screenData_normalized)
"screenData_normalized"
