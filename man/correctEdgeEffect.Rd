% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correctEdgeEffect.R
\name{correctEdgeEffect}
\alias{correctEdgeEffect}
\title{Estimate and correct edge effect for each screen plate}
\usage{
correctEdgeEffect(
  screenData,
  estimateMethod = "loess",
  useNeg = TRUE,
  useLowConcentrations = 0,
  span = 1,
  correctMethod = "bliss",
  exclude = c()
)
}
\arguments{
\item{screenData}{a data frame containing screening results generated by \code{readScreen()} function. In order for the edge effect estimate algorithm to work,
plate normalization needs to be performed first.}

\item{estimateMethod}{a character string specifying the method used for estimating the shape and strength of edge effect.
Currently two methods are supported: 'loess', which uses a local regression model to estimate edge effect;
'sigmoid', an experimental feature, which uses use a 2D-sigmoid model to estimate edge effect.
The "loess" method is suitable for both regular incubation effect, such as edge effect, and irregular incubation effect. But sometimes it tends to over-smooth the edge effect
when the edge effect is rectangular. In this case, "sigmoid" method can better capture the shape of edge effect. But the "sigmoid" method can not deal with incubation effect with irregular shape.}

\item{useNeg}{a logical value. If TRUE, only negative controls are used to estimate edge effect. Otherwise, all the wells are used for edge effect estimation.}

\item{useLowConcentrations}{a integer value. In addition to negative controls,
sample wells with low concentrations of compounds can be considered as 'de facto' negative controls and included in edge effect estimation.
The default value is 0, means no sample wells are used. If the value n >= 1, the n lowest concentrations will be used.}

\item{span}{a numeric value. The span parameter for \code{loess} method, which controls the degree of smoothing.}

\item{correctMethod}{a character string specifying the method used for subtracting the edge effect. Currently two methods 'bliss' and 'linear' are supported.
If \code{bliss} is selected, the corrected values are simply the normalized signal intensity values divided by the estimated incubation/edge effect intensity.
If \code{linear} is selected, the correction process will be smoothed by a linear function depends on the original normalized signal intensity. This option is
suitable for the situation when the edge effect is very strong and results in instability in the correction process.}

\item{exclude}{a list of samples that should be excluded when performing edge effect correction}
}
\value{
This function will add one column, \code{normVal.cor}, which is the normalized viability corrected for incubation/edge effect, to the input data frame.
}
\description{
This function subtract the incubation (or edge) effect from the observed signal intensity to recover the true signal intensity.
If the incubation effect has already been estimated by \code{fitEdgeEffect()} function, the values in the edgeFactor column will be used.
Otherwise, this function will first call \code{fitEdgeEffect()} function to estimate the edge effect.
}
\examples{
# load drug screen dataset
data('screenData_normalized')
screenData_corrected <- correctEdgeEffect(screenData_normalized, estimateMethod = 'loess',
    useNeg = TRUE, useLowConcentrations = 0, span =1, correctMethod = 'bliss')
}
