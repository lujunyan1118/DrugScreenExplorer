#' Create a template for well annotation file
#'
#' This function creates a tab-separated values (tsv) file required for annotating wells on drug screening plate, for example, the name of the drug or drug concentration in each well.
#' Users can use the file created by this function as a starting point for preparing the well annotation file required by the \code{readScreen()} function.
#' Users can also prepare the well annotation file from scratch, as long as it follows the format requirement. More information is available in the package vignette.
#'
#' @param file the output file name and directory for the sample annotation file
#' @param colNum the number of columns on the plate, will be labeled by numeric numbers
#' @param rowNum the number of rows on the plate, will be labeled by alphabets
#' @param entries  a vector of characters to use as the column names in the drug file (optional)
#' @param platePerSample a numeric value indicates the number of plates used for each sample in the screen. The default value is 1.
#' @return This function creates a tab-separated values (tsv) file in the working directory.
#' @import tibble dplyr
#' @export
#' @examples
#' # create a template for well annotation file
#' createWellInput(file = 'wellAnno.tsv', colNum = 24,
#'     rowNum = 16, entries = c('name', 'concentration'))
#' # see vignette for more instructions


createWellInput <- function(file, colNum, rowNum, entries = c("name", "concentration"), platePerSample = 1) {
  rowID <- genRowIDs(rowNum)
  colID <- genColIDs(colNum)
  rowID <- rep(rowID, each = colNum)
  colID <- rep(colID, times = rowNum)
  outTable <- tibble::tibble(wellID = paste0(rowID, colID))

  for (eachEntry in entries) {
    outTable[[eachEntry]] <- ""
  }

  if (platePerSample > 1) {
    newTab <- lapply(seq(platePerSample), function(x) {
      dplyr::mutate(outTable, plateID = x)
    }) %>%
      dplyr::bind_rows()
    outTable <- newTab
  }

  write_delim(outTable, file = file, delim = "\t")

}
