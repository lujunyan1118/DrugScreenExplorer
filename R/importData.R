#Collection of functions for import screen data and annotations

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Create plate annotation file
#'
#' This function creates a csv file required for annotating plates, based on the raw data files
#'
#' @param rawDir A character string specifying the directory of the raw data folder
#' @param file The output file name and path for the plate annotation file
#' @param entries A vector of characters you wish to use as the column names in the plate annotation file (optional)
#' @param batchAsFolder A logical value. Whether the raw data from different batches are stored in different subfolders. Default is FALSE.
#' @param csvFormat A character sting indicating the format of the csv file, either "csv" (comma separated) or "csv2" (semicolon separated). Default is "csv2".
#' @import tidyverse tools
#' @export

createPlateInput <- function(rawDir, file, entries = c(), batchAsFolder = FALSE, csvFormat = "csv2") {
  if (!batchAsFolder) {
    allFiles <- list.files(path = rawDir, recursive = TRUE)
    outTable <- tibble(fileName = basename(allFiles))
  } else {
    #raw data from different batches are stored in subfolders
    batches <- list.dirs(rawDir,full.names = FALSE)
    batches <- sort(batches[batches != ""])
    outTable <- lapply(seq(length(batches)), function(batch) {
      allFiles <- list.files(path = file.path(rawDir, batches[batch]))
      tibble(fileName = allFiles, batch = batch)
    }) %>% bind_rows()

  }

  #check for duplicated file names
  if (sum(duplicated(outTable$fileName)) > 0 ) stop("Duplicated file names found")

  for (eachEntry in entries) {
    outTable <- add_column(outTable, !!eachEntry := "")
  }

  if (csvFormat == "csv") {
    write.csv(outTable, file = file, row.names = FALSE)
  } else {
    write.csv2(outTable, file = file, row.names = FALSE)
  }
}

#' Create well annotation file
#'
#' This function creates a csv file required for annotating wells on drug screening plate
#'
#' @param file The output file name and directory for the sample annotation file
#' @param colNum The number of columns on the plate, will be labelled by number
#' @param rowNum The number of rows on the plate, will be labelled by alphabets
#' @param entries  A vector of characters you wish to use as the column names in the drug file (optional)
#' @param platePerSample The number of plates used for each sample in the screen. Default is 1.
#' @import tidyverse
#' @export

createWellInput <- function(file, colNum, rowNum, entries = c("name","concentration"), platePerSample = 1, csvFormat = "csv2") {
  rowID <- genRowIDs(rowNum)
  colID <- genColIDs(colNum)
  rowID <- rep(rowID, each = colNum)
  colID <- rep(colID, times = rowNum)
  outTable <- tibble(wellID = paste0(rowID, colID))

  for (eachEntry in entries) {
    outTable <- add_column(outTable, !!eachEntry := "")
  }

  if (platePerSample > 1) {
    newTab <- lapply(seq(platePerSample), function(x) {
      mutate(outTable, plateID =x)
    }) %>% bind_rows()
    outTable <- newTab
  }

  if (csvFormat == "csv") {
    write.csv(outTable, file = file, row.names = FALSE)
  } else {
    write.csv2(outTable, file = file, row.names = FALSE)
  }

}

#' Read the screening result from one plate
#'
#' This function creates a tidy table containing the screening result from one plate, mainly used by the readScreen function
#'
#' @param plateFile A characther string specifying the file name for the plate result
#' @param rowRange Either An integer vector with length of two or an integer value. If it's an vector, the range of the lines that contain the result values are speficified. Otherwise, only the starting line is indicated.
#' @param colRange Either An integer vector with length of two or an integer value. If it's an vector, the range of the columns (separated by delimiters) that contain the result values are speficified. Otherwise, only the starting column is indicated.
#' @param sep A character string specifying the delimiters for extracting values in the text file. A regular expression can be used. Default value is using both tab and semicolon as delimiter.
#' @param commaAsDecimal A logical values, whether commas is used as decimal in the screen value. Default if FALSE
#' @import tidyverse
#' @export
#'

readPlate <- function(plateFile,rowRange = c(3,18), colRange = 2, sep = "[;\t]",
                      commaAsDecimal = FALSE) {
  txt <-readLines(plateFile)

  if(length(rowRange) == 1) {
    rowRange <- c(rowRange, length(txt))
  }

  sp <- lapply(seq(rowRange[1],rowRange[2]),  function(x){
    if (commaAsDecimal)
      txtSp <- gsub("[,]",".",txt[x]) else
        txtSp <- txt[x] #only "," or "." as decimal
    txtSp <- strsplit(txtSp, split = sep)[[1]]
    txtSp <- txtSp[txtSp != ""]
    if (length(colRange) == 1) colRange <- c(colRange, length(txtSp))
    txtSp <- txtSp[colRange[1]:colRange[2]]
  })

  #check whether all the rows have the same length.
  rowNum <- length(sp)
  colNum <- unique(sapply(sp, length))
  if(length(colNum) != 1) {
    stop("Numbers of the values in each row are different.")
  }

  rowID <- genRowIDs(rowNum)
  colID <- genColIDs(colNum)
  rowID <- rep(rowID, each = colNum)
  colID <- rep(colID, times = rowNum)

  df <- tibble(wellID = paste0(rowID, colID),rowID = rowID,
                               colID = colID, value = as.numeric(unlist(sp)))
  stopifnot(!any(is.na(df$value)))

  return(df)
}

#' Read the screening results from a whole experiment
#'
#' This function creates a tidy table containing the screening result from a whole experiment. Drugs and samples are annotated based on drug annotation and sample annotation file
#'
#' @param rawDir A character string specifying the directory of the raw data folder.
#' @param plateAnnotationFile A character string specifying the plate annotation file.
#' @param wellAnnotationFile A character string specifying the well annotation file.
#' @param rowRange Either An integer vector with length of two or an integer value. If it's an vector, the range of the lines that contain the result values are speficified. Otherwise, only the starting line is indicated.
#' @param colRange Either An integer vector with length of two or an integer value. If it's an vector, the range of the columns (separated by delimiters) that contain the result values are speficified. Otherwise, only the starting column is indicated.
#' @param sep A character string specifying the delimiters for extracting values in the text file. A regular expression can be used. Default value is using both tab and semicolon as delimiter.
#' @param negWell A character vector or a character string of well IDs or drug names (if the "name" column is present in well annotation), which specifies the negative control wells
#' @param posWell A character vector or a character string of well IDs or drug names (if the "name" column is present in well annotation), which specifies the positive control wells
#' @param normalization A logical value, whether to calculate percent inhibition based on control wells on the plate. Default is FALSE
#' @param method A character string, either "negatives" or "both", specifying the method for calculating percent inhibition values. If "negatives" is used, each measurement on the plate will be divided by the median of the values of negative control wells on the same plate. If "both" is used, a "normalized percent inhbition (NPI)" is applied in a per-plate basis. For an inhibition assay, this method divides the difference between each measurement on a plate and the median measurement of the positive controls on that plate by the difference between the median mesurement of negative controls and positive controls on that well. The default value is "negatives".
#' @param discardLayer An integer values. To avoid the impact of incubation effect on the negative controls, a certain number of the exterior on the plate can be ignored when summarising the measurement of controls. 0 means do not discard any control wells.
#' @param commaAsDecimal A logical value, whether commas is used as decimal in the screen value. Default if FALSE.
#' @param csvFormat A character string specifying the csv file type, either "csv" or "csv2". Default is csv2
#' @import tidyverse
#' @export
#'
readScreen <- function(rawDir, plateAnnotationFile, wellAnnotationFile,
                       rowRange, colRange, sep,
                       negWell = c(),posWell = c(),
                       normalization = FALSE, method = "negatives", discardLayer = 0,
                       commaAsDecimal = FALSE, csvFormat = "csv2") {

  rawFileList <- list.files(rawDir, full.names = TRUE, recursive = TRUE)

  #check if there are duplicated file names
  allFiles <- sapply(rawFileList, basename)
  if (sum(duplicated(allFiles)) >0 ) stop("Duplicated file names found")

  #read screening data plate by plate
  screenData <- lapply(rawFileList, function(fileName) {

    plateData <- tryCatch({
      readPlate(plateFile = fileName, rowRange, colRange, sep, commaAsDecimal) %>%
      mutate(fileName = basename(fileName)) },
      error = function(e) {
        stop(sprintf("Error encountered when reading %s", fileName))
      })

  }) %>% bind_rows()

  #read sample annotation files
  if (csvFormat == "csv2")
    sampleAnno <- read_csv2(plateAnnotationFile) else
      sampleAnno <- read_csv(plateAnnotationFile)

  #add sample annotations to screenData
  screenData <- left_join(screenData,  sampleAnno, by = "fileName")

  #remove extensions in the file names
  screenData <- mutate(screenData, fileName = file_path_sans_ext(fileName))

  #read drug annotation files
  if (csvFormat == "csv2")
    drugAnno <- read_csv2(wellAnnotationFile) else
       drugAnno <- read_csv(wellAnnotationFile)

  #add drug annotation information
  if ("plateID" %in% colnames(drugAnno))
    screenData <- left_join(screenData, drugAnno, by = c("wellID", "plateID"))
  else  screenData <- left_join(screenData, drugAnno, by = "wellID")

  #add well type information
  if (length(intersect(negWell, posWell)) != 0) stop("Overlap in postive and negative control well list")
  if ("name" %in% colnames(screenData)) {
    screenData <- mutate(screenData, wellType = ifelse(wellID %in% negWell | name %in% negWell, "neg","sample")) %>%
      mutate(wellType = ifelse(wellID %in% posWell | name %in% posWell, "pos", wellType))
  } else {
    screenData <- mutate(screenData, wellType = ifelse(wellID %in% negWell, "neg","sample")) %>%
      mutate(wellType = ifelse(wellID %in% posWell, "pos", wellType))
  }

  #normalizing plates
  if (normalization) {
    screenData <- normalizePlate(screenData = screenData, method = method, discardLayer = discardLayer)
  }

  #final adjustment
  ## batch as factor
  if ("batch" %in% colnames(screenData)) screenData$batch <- as.factor(screenData$batch)
    screenData <- mutate_if(screenData, is.character, as.factor) %>% ungroup(screenData)
  return(screenData)
}

#' Per-plate normalization
#'
#' This function calculates the percent inhibition values for each plate based on the control wells on the same plate.
#' @param screenData A tidytable that contains the output from a readScreen function.
#' @param method A character string, either "negatives" or "both", specifying the method for calculating percent inhibition values. If "negatives" is used, each measurement on the plate will be divided by the median of the values of negative control wells on the same plate. If "both" is used, a "normalized percent inhbition (NPI)" is applied in a per-plate basis. For an inhibition assay, this method divides the difference between each measurement on a plate and the median measurement of the positive controls on that plate by the difference between the median mesurement of negative controls and positive controls on that well. The default value is "negatives".
#' @param discardLayer An integer values. To avoid the impact of incubation effect on the negative controls, a certain number of the exterior on the plate can be ignored when summarising the measurement of controls. 0 means do not discard any control wells.
#' @import tidyverse
#' @export

normalizePlate <- function(screenData, method = "negatives", discardLayer = 0) {

  if (discardLayer > 0) {
    #generate a list of wellIDs that need to be ignored when summarising control measurements
    nCol <- length(unique(screenData$colID))
    nRow <- length(unique(screenData$rowID))
    disCol <- genColIDs(nCol)[c(seq(1,discardLayer), seq(nCol-discardLayer+1, nCol))]
    disRow <- genRowIDs(nRow)[c(seq(1,discardLayer), seq(nRow-discardLayer+1, nRow))]
  } else {
    disCol <- c()
    disRow <- c()
  }

  normOnePlate <- function(onePlate) {
    negMed <- median(filter(onePlate, wellType == "neg", ! colID %in% disCol, !rowID %in% disRow)$value, na.rm = TRUE)
    posMed <- median(filter(onePlate, wellType == "pos", ! colID %in% disCol, !rowID %in% disRow )$value, na.rm = TRUE)
    if (method == "negatives") {
      if (is.na(negMed)) {
        stop("No negative control wells found on the speficied region of the plate")
      } else {
        onePlate <- mutate(onePlate, normVal = value/negMed)
      }
    } else if (method == "NPI") {
      if (is.na(negMed) | is.na(posMed)) {
        stop("No negative and poitive control wells found on the specified region of the plate")
      } else {
        onePlate <- mutate(onePlate, normVal = (value - posMed)/(negMed - posMed))
      }
    } else stop('Please choose a normalization method, either "negatives" or "NPI"')
    return(onePlate)
  }

  group_by(screenData, fileName) %>% do(normOnePlate(.)) %>% ungroup()
}
