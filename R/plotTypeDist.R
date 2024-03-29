#' Plotting raw signal distribution for each well type
#'
#' A function to generate box plots to show the raw signal distribution for each well type on all or selected plates.
#' For this function to work, the "wellType" column must be present in the data frame.
#'
#' @param screenData a data frame containing the screen data generated by readScreen() function
#' @param plotPlate a character string or a vector of character strings specifying the file name of the plates to be included in the plots.
#' If all plates should be included, 'all' can be used. The default value is 'all'.
#' @param ifLog10 a logical value, whether to perform a log10 transformation on the values.
#' @param pdfName a character string specifying the output pdf file name.
#' If not specified, a list containing the ggplot objects will be returned. Default is NULL.
#' @export
#' @import ggplot2 dplyr
#' @return If the argument \code{pdfName = NULL}, this function returns a list of ggplot objects. If a file name is specified \
#' for the \code{pdfName} argument, a pdf file that contains the plots will be created in the working directory.
#' @examples
#' # load processed data
#' data('screenData')
#'
#' # plot raw ATP count distribution for each well type
#' plotTypeDist(screenData)
#' # Please see the vignette for more information.

plotTypeDist <- function(screenData, plotPlate = "all", ifLog10 = FALSE, pdfName = NULL) {

  if (!"wellType" %in% colnames(screenData))
    stop("No well type annotation found")

  if (plotPlate == "all") {
    plateList <- unique(screenData$fileName)

  } else {
    plateList <- intersect(unique(screenData$fileName), plotPlate)

  }
  if (length(plateList) == 0)
    stop("No plate found")

  pList <- lapply(plateList, function(plateName) {
    plotTab <- dplyr::filter(screenData, fileName == plateName) %>%
      dplyr::mutate(wellType = factor(wellType))
    if (length(unique(plotTab$wellType)) <= 1)
      stop("The well types on plate is less than 1, nothing to plot")

    if (ifLog10) {
      plotTab <- dplyr::mutate(plotTab, value = log10(value))
      yLabText <- ylab("log10(raw counts)")

    } else {
      yLabText <- ylab("raw counts")

    }

    ggplot(plotTab, aes(x = wellType, y = value, fill = wellType)) + geom_point() +
      geom_boxplot(alpha = 0.5) + theme_bw() +
      theme(legend.position = "none",plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
            axis.text.x = element_text(face = "bold",size = 10),
            axis.text.y = element_text(size = 10), plot.margin = margin(5,5, 5, 5)) +
      yLabText + xlab("well types") + ggtitle(plateName)
  })
  names(pList) <- plateList

  if (is.null(pdfName)) {
    return(pList)

  } else {
    makepdf(pList, pdfName, ncol = 2, nrow = 2, 8, 8)

  }
}
