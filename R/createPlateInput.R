#' Create a template for plate annotation file
#'
#' This function creates a tab-separated values (tsv) file required for annotating plates, for example, the sample ID or patient ID for each plate.
#' Users can use the file created by this function as a starting point for preparing the plate annotation file required by the \code{readScreen()} function.
#' Users can also prepare the plate annotation file from scratch, as long as it follows the format requirement. More information is available in the package vignette.
#'
#' @param rawDir a character string specifying the directory of the raw data folder
#' @param file the output file name and path for the plate annotation file
#' @param entries a vector of character strings to use as the column names in the plate annotation file (optional)
#' @param batchAsFolder a logical value, whether the raw data from different batches are stored in different sub-folders.
#' The default value is \code{FALSE}. If \code{TRUE}, the folder names will be used as batch identifiers automatically.
#' @return This function creates a tab-separated values (tsv) file in the working directory.
#' @import tibble dplyr readr
#' @export
#' @examples
#' # create a template for plate annotation file
#' rawFolder <- system.file('testData/rawData', package = 'DrugScreenExplorer')
#' createPlateInput(rawDir = rawFolder, file = 'plateAnno.tsv',
#'     entries = c('sampleID', 'patientID'))
#' # see vignette for more instructions

createPlateInput <- function(rawDir, file, entries = c(), batchAsFolder = FALSE) {
  if (!batchAsFolder) {
    allFiles <- list.files(path = rawDir, recursive = TRUE)
    outTable <- tibble::tibble(fileName = basename(allFiles))
  } else {
    # raw data from different batches are stored in subfolders
    batches <- list.dirs(rawDir, full.names = FALSE)
    batches <- sort(batches[batches != ""])

    outTable <- lapply(seq(length(batches)), function(batch) {
      allFiles <- list.files(path = file.path(rawDir, batches[batch]))
      tibble::tibble(fileName = allFiles, batch = batch)
    }) %>%
      dplyr::bind_rows()

  }

  # check for duplicated file names
  if (sum(duplicated(outTable$fileName)) > 0)
    stop("Duplicated file names found")

  for (eachEntry in entries) {
    outTable[[eachEntry]] <- ""
  }

  readr::write_delim(outTable, file = file, delim = "\t")

}
