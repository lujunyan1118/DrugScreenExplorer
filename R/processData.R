# Collection of functions for processing of screen data

#' Model object for fitting the IC50 curve
#'
#' Use the drc package to perform IC50 fit. Can be directly used for geom_smooth() in ggplot2
#' @param formula Formula for the curve fitting.
#' @param data A data frame contain the raw concentration and the viability value. The viability should not be the percent viability value.
#' @param weigths Not used, mainly for geom_smooth() purpose
#' @param ... Parameters passed to logLogisticRegression()
#' @export
#' @import drc
fitIC50 <- function(formula, data = NULL, weights, ...) {
  if (! is.null(data) ) {
    modelFrame <- model.frame(formula, data)
  } else {
    modelFrame <- model.frame(formula)
  }
  parm_fit <- drm(modelFrame, fct = LL2.3u(), ...)
  newModel <- list(model = modelFrame, formula = formula, parm_fit = parm_fit)
  class(newModel) <- "fitIC50"
  return(newModel)
}


#' Predicted values based on IC50 fit
#'
#' Generic function for ic50 class generated from fitIC50 function.
#' @param object Object of class inheriting from "fitIC50"
#' @param newdata An optional data frame in which to look for variables with which to predict.If omitted, the fitted values are used.
#' @param se.fit Not used, mainly for geom_smooth purpose
#' @param level Not used, mainly for geom_smooth purpose
#' @param interval Not used, mainly for geom_smooth purpose
#' @export
#'
predict.fitIC50 <- function(object, newdata = NULL, se.fit = FALSE, level = 0.95 , interval = c("none", "confidence", "prediction")) {

  if (is.null(newdata))
    newdata <- object$model else
      newdata <- newdata

    parm_fit <- object$parm_fit

    res <- predict(parm_fit, newdata)

    return(res)
}
