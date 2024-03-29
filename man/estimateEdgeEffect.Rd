% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateEdgeEffect.R
\name{estimateEdgeEffect}
\alias{estimateEdgeEffect}
\title{Estimate edge effect intensity}
\usage{
estimateEdgeEffect(
  screenData,
  nLayer = 1,
  onlyNeg = TRUE,
  identifier = "fileName",
  pdfName = NULL
)
}
\arguments{
\item{screenData}{a data frame containing the screening results generated by \code{readScreen()} function}

\item{nLayer}{an integer value (larger than 1), which indicates up to how many layers on the edge of the screen are considered as edges. The default value is 1.}

\item{onlyNeg}{a logical value, whether only wells annotated as negative controls are used for edge effect estimation.
If \code{FALSE}, all the measurements on each plate are used for edge effect estimation. The default value is \code{TRUE}.}

\item{identifier}{a character string specifying the column in the input data frame that should be used as identifiers (x-axis) in the plot.
Default value is "fileName". If one identifier matches several plates, error bars indicating standard deviations will be shown in the plot.}

\item{pdfName}{either NULL or a character string specifying the output file name for the plot in pdf format. If NULL, no plot will be created.}
}
\value{
Depends on the number \emph{N} specified with the \code{nLayer} argument, this function will add \emph{N} columns with the names edgeRatio1, edgeRatio2, ... edgeRation\emph{N}. Each column contains the estimated edge effect intensity when considering the \emph{N} outer layers on the plate as edge layers.
If a file name is specified for the \code{pdfName} argument, a pdf file that contains a bar plot visualizing the edge effect intensity will be created in the working directly.
}
\description{
This function estimates the intensity of edge effect on each screen plate by comparing the median signal intensities on the edge to those on the inner plate.
The results will be shown as a bar plot with plate IDs as x-axis and edge effect intensities as y-axis.
}
\examples{
# load processed data
data('screenData_normalized')

# estimate the intensity of edge effect.
# The wells from the first outside layer on the plate will be defined as the edge wells.
estimateEdgeEffect(screenData_normalized, nLayer =1)

# Please see the vignette for more information.
}
