% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importData.R
\name{readPlate}
\alias{readPlate}
\title{Read the screening result from one plate}
\usage{
readPlate(plateFile, rowRange = c(3, 18), colRange = 2,
  sep = "[;\\t]", commaAsDecimal = FALSE)
}
\arguments{
\item{plateFile}{A characther string specifying the file name for the plate result}

\item{rowRange}{Either An integer vector with length of two or an integer value. If it's an vector, the range of the lines that contain the result values are speficified. Otherwise, only the starting line is indicated.}

\item{colRange}{Either An integer vector with length of two or an integer value. If it's an vector, the range of the columns (separated by delimiters) that contain the result values are speficified. Otherwise, only the starting column is indicated.}

\item{sep}{A character string specifying the delimiters for extracting values in the text file. A regular expression can be used. Default value is using both tab and semicolon as delimiter.}

\item{commaAsDecimal}{A logical values, whether commas is used as decimal in the screen value. Default if FALSE}
}
\description{
This function creates a tidy table containing the screening result from one plate, mainly used by the readScreen function
}
