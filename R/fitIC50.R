# Fitting 4 parameter log-logistic (4PL) models.
#'
#' Use the dr4pl package to perform IC50 curve fitting. Can be directly used for geom_smooth() in ggplot2.
#'
#' @param formula formula for the curve fitting.
#' @param data a data frame that contains the  concentration values and the normalized viability values or percent inhibition values.
#' @param weights not used, mainly for geom_smooth() purpose
#' @param logDose a numeric value specifying the base of log transformation for concentration,
#' if the concentration is log transformed. The default value is NULL (not log transformed).
#' @param ... further arguments to be passed to dr4pl function
#' @export
#' @import dr4pl
#' @return fitIC50 returns an object of class 'fitIC50', which is a list containing the following components:
#' \item{model}{the input dose-response table}
#' \item{formula}{the formula used for the model fitting}
#' \item{parm_fit}{the fitted parameters for a 4-paramater logistic model}
#' \item{logDose}{a numeric value specifying the base of log transformation for concentration,
#' if the concentration is log transformed. The default value is NULL (not log transformed).}
#' @examples
#' # create an example dose-response table
#' doseTab <- data.frame(viability = c(1.0, 0.95, 0.5, 0.2, 0.1),
#'                       concentration = c(0.001, 0.01, 0.1, 1, 10))
#'
#' # fit a four-parameter logistic model using dr4pl package
#' fitIC50(viability ~ concentration, data = doseTab)
#'
fitIC50 <- function(formula, data = NULL, weights = NULL, logDose = NULL, ...) {
  if (!is.null(data)) {
    modelFrame <- model.frame(formula, data)
  } else {
    modelFrame <- model.frame(formula)
  }

  if (!is.null(logDose)) {
    modelFrame[, 2] <- logDose^modelFrame[, 2]
  }

  parm_fit <- dr4pl(modelFrame[, 2], modelFrame[, 1], ...)
  newModel <- list(model = modelFrame, formula = formula, parm_fit = parm_fit, logDose = logDose)

  class(newModel) <- "fitIC50"
  return(newModel)
}
