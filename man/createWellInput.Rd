% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createWellInput.R
\name{createWellInput}
\alias{createWellInput}
\title{Create a template for well annotation file}
\usage{
createWellInput(
  file,
  colNum,
  rowNum,
  entries = c("name", "concentration"),
  platePerSample = 1
)
}
\arguments{
\item{file}{the output file name and directory for the sample annotation file}

\item{colNum}{the number of columns on the plate, will be labeled by numeric numbers}

\item{rowNum}{the number of rows on the plate, will be labeled by alphabets}

\item{entries}{a vector of characters to use as the column names in the drug file (optional)}

\item{platePerSample}{a numeric value indicates the number of plates used for each sample in the screen. The default value is 1.}
}
\value{
This function creates a tab-separated values (tsv) file in the working directory.
}
\description{
This function creates a tab-separated values (tsv) file required for annotating wells on drug screening plate, for example, the name of the drug or drug concentration in each well.
Users can use the file created by this function as a starting point for preparing the well annotation file required by the \code{readScreen()} function.
Users can also prepare the well annotation file from scratch, as long as it follows the format requirement. More information is available in the package vignette.
}
\examples{
# create a template for well annotation file
createWellInput(file = 'wellAnno.tsv', colNum = 24,
    rowNum = 16, entries = c('name', 'concentration'))
# see vignette for more instructions
}
