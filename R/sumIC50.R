#' Summarise drug effect across concentrations using a 4-parameter log-logistic model
#'
#' This function will fit a 4-parameter log-logistic model using dr4pl package based on the input dose-response table and return the model parameters as well as integrated area under the fitted sigmoid curve.
#' @param formula formula for the curve fitting.
#' @param data a data frame containing the raw concentration and the viability value.
#' @param minConc lower bound of concentrations for the integration of area under curve, default value is 0.
#' @param maxConc upper bound of concentrations for the integration of area under curve, default value is 15.
#' @param n number of bins for integration.
#' @param ... further arguments to be passed to dr4pl function.
#' @export
#' @import dr4pl
#' @return a dataframe with only one row and four columns:
#' \item{UpperLimit}{the upper limit of the fitted logistic model}
#' \item{IC50}{the half maximal inhibitory concentration (IC50),
#' or the half maximal effective concentration (EC50), depends on the type of input data}
#' \item{Slope}{the slope of the fitted logistic model}
#' \item{lowerLimit}{the upper limit of the fitted logistic model}
#' \item{AUC}{the area under the fitted sigmoid doese-response curve}
#' If the model fitting is failed, all values in the dataframe will be NA.
#'
#' @examples
#' # create an example dose-response table
#' doseTab <- data.frame(viability = c(1.0, 0.95, 0.5, 0.2, 0.1),
#'                       concentration = c(0.001, 0.01, 0.1, 1, 10))
#'
#' # get the parameters
#' sumIC50(viability ~ concentration, data = doseTab)

sumIC50 <- function(formula, data = NULL, minConc = 0, maxConc = 15, n = 100, ...) {
  res <- tryCatch({
    # fit sigmoid model
    mm <- fitIC50(formula, data, ...)

    # integrate AUC
    concSeq <- data.frame(conc = seq(minConc, maxConc, length.out = n))
    valueConc <- data.frame(viab = predict.fitIC50(mm, concSeq), conc = concSeq[, 1])
    aucVal <- calcAUC(valueConc$viab, valueConc$conc)

    res <- data.frame(t(mm$parm_fit$parameters))
    colnames(res) <- c("UpperLimit", "IC50", "Slope", "LowerLimit")
    res$AUC <- aucVal
    res
  }, error = function(err) {
    warning("Curve fitting failed, NA values generated. ")
    data.frame(t(structure(rep(NA, 5), names = c("UpperLimit", "IC50", "Slope", "LowerLimit", "AUC"))))
  })

  return(res)
}
