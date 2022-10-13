#' Read and assemble a drug screen dataset
#'
#' This function creates a data frame that contains the screening result from a whole drug screen experiment.
#' Drugs and samples are annotated based on the drug annotation and sample annotation files.
#'
#' @param rawDir a character string specifying the location of the raw data folder.
#' @param plateAnnotationFile a character string specifying the location of the plate annotation file.
#' @param wellAnnotationFile a character string specifying the location of the well annotation file.
#' @param rowRange either an integer vector with length of two or an integer value.
#' If it's an vector, it should specify the range of the lines that contain the actual values.
#' If it's an integer values, only the starting line is indicated and the end line is the end of the file.
#' @param colRange either an integer vector with length of two or an integer value.
#' If it's an vector, it should specify the range of the columns (separated by delimiters) that contain the actual values.
#' If it's an integer, only the starting column is indicated and the end column is the last column in the file.
#' @param sep a character string specifying the delimiters for extracting values in the text file.
#' A regular expression can be used. The default value is using either tab or semicolon as delimiter.
#' @param negWell a character vector or a character string of well IDs or drug names
#' (if the 'name' column is present in well annotation), which specifies the negative control wells
#' @param posWell a character vector or a character string of well IDs or drug names
#' (if the 'name' column is present in well annotation), which specifies the positive control wells
#' @param normalization a logical value to indicate whether to calculate percent inhibition
#' based on control wells on the plate. Default is FALSE.
#' @param method a character string, either 'negatives' or 'both',
#' specifying the method for calculating percent inhibition values. If 'negatives' is used,
#' each measurement on the plate will be divided by the median of the values of negative control wells on the same plate.
#' If 'both' is used, a 'normalized percent inhibition (NPI)' is applied in a per-plate basis.
#' For an inhibition assay, this method divides the difference between each measurement on a plate and
#' the median measurement of the positive controls on that plate by the difference between the median measurement of
#' negative controls and positive controls on that well. The default value is 'negatives'.
#' @param discardLayer an integer values. To avoid the impact of incubation effect on the negative controls,
#' a certain number of the exterior wells on the plate can be ignored when summarising the measurement of controls.
#' 0 means do not discard any control wells. A value above zero indicates the \code{N} layers of wells on the edge should be discarded.
#' @param commaAsDecimal a logical value, whether commas is used as decimal in the input data. Default if FALSE.
#' @param allowNA a logical value, whether to allow NA (missing or empty) values in the input data.
#' @param bySample a logic value, indicates whether to normalize the viability by sample rather than by plate.
#' This is for situations where multiple samples are on the same plate.
#' @import tibble dplyr readr stringr tools
#' @export
#' @return A data frame with each row represent a well on a plate and the following columns:
#' \item{wellID}{well identifier}
#' \item{rowID}{row identifier}
#' \item{colID}{column identifier}
#' \item{value}{raw signal intensity from a plate reader}
#' \item{fileName}{input files names (without extension names) for the plates}
#' \item{name}{drug name, if this column is present in the well annotation file}
#' \item{concentration}{drug concentration, if this column is present in the well annotation file}
#' \item{wellType}{the type of wells, can be either neg (negative controls), pos (positive controls) or sample (drug wells)}
#' \item{normVal}{normalized viability based on control wells on the plate, if the normalization step is performed by setting \code{normalization = TRUE}}
#' Depends on the user-specified plate and well annotations, other columns, such as sampleID, patientID,
#' batch and so on, in the annotation files can also be present in the output table.
#' @examples
#' rawFolder <- system.file('testData/rawData', package = 'DrugScreenExplorer')
#' plateFile <- system.file('testData/plateAnno/plateAnno.tsv', package = 'DrugScreenExplorer')
#' wellFile <- system.file('testData/wellAnno/wellAnno.tsv', package = 'DrugScreenExplorer')
#' screenData <- readScreen(rawDir = rawFolder, plateAnnotationFile = plateFile,
#'                         wellAnnotationFile = wellFile,
#'                         rowRange = c(3,18), colRange = 2,
#'                         sep = '[;\t]', commaAsDecimal = FALSE,
#'                         negWell <- c('DMSO','PBS'), posWell = c())
#' # see vignette for more instructions


readScreen <- function(rawDir, plateAnnotationFile, wellAnnotationFile, rowRange, colRange, sep,
                       negWell = c(), posWell = c(), normalization = FALSE, method = "negatives", discardLayer = 0,
                       commaAsDecimal = FALSE, allowNA = FALSE, bySample = FALSE) {

  rawFileList <- list.files(rawDir, full.names = TRUE, recursive = TRUE)

  # read sample annotation files
  if (commaAsDecimal) {
    deciMark <- ","
  } else {
    deciMark <- "."
  }

  sampleAnno <- read_delim(plateAnnotationFile, delim = "\t", locale = locale(decimal_mark = deciMark), col_types = cols())

  # check if there are duplicated file names
  allFiles <- vapply(rawFileList, basename, "c")

  if (sum(duplicated(allFiles)) > 0) stop("Duplicated file names found")

  # read screening data plate by plate
  screenData <- lapply(rawFileList, function(fileName) {

    plateData <- tryCatch({
      readPlate(plateFile = fileName, rowRange = rowRange, colRange = colRange, sep = sep,
                commaAsDecimal = commaAsDecimal, allowNA = allowNA) %>%
        dplyr::mutate(fileName = basename(fileName))
    }, error = function(e) {
      stop(sprintf("Error encountered when reading %s", fileName))
    })

  }) %>% bind_rows()


  # add sample annotations to screenData
  screenData <- left_join(screenData, sampleAnno, by = "fileName")

  # remove extensions in the file names
  screenData <- dplyr::mutate(screenData, fileName = tools::file_path_sans_ext(fileName))

  # read drug annotation files
  drugAnno <- read_delim(wellAnnotationFile, delim = "\t", locale = locale(decimal_mark = deciMark), col_types = cols())

  # add drug annotation information
  if ("plateID" %in% colnames(drugAnno)) {
    screenData <- left_join(screenData, drugAnno, by = c("wellID", "plateID"))
  } else {
    screenData <- left_join(screenData, drugAnno, by = "wellID")
  }

  # add well type information
  if (length(intersect(negWell, posWell)) != 0)
    stop("Overlap in postive and negative control well list")

  if ("name" %in% colnames(screenData)) {
    screenData <- dplyr::mutate(screenData, wellType = ifelse(wellID %in% negWell | name %in%
                                                                negWell, "neg", "sample")) %>%
      dplyr::mutate(wellType = ifelse(wellID %in% posWell | name %in% posWell, "pos", wellType))
  } else {
    screenData <- dplyr::mutate(screenData, wellType = ifelse(wellID %in% negWell, "neg", "sample")) %>%
      dplyr::mutate(wellType = ifelse(wellID %in% posWell, "pos", wellType))
  }

  if (normalization) { # normalizing plates
    screenData <- normalizePlate(screenData = screenData, method = method,
                                 discardLayer = discardLayer,bySample = bySample)
  }

  if ("batch" %in% colnames(screenData)) {     # batch as factor
    screenData$batch <- factor(screenData$batch)
  }

  if ("concentration" %in% colnames(screenData)) {     # concentration as numeric
    screenData$concentration <- as.numeric(as.character(screenData$concentration))
  }

  screenData <- ungroup(screenData)  #%>% mutate_if(is.character, factor)

  return(screenData)
}
