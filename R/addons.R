# Some helper functions

## Function to make pdf for a list of ggplot objects
makepdf <- function(x, name, ncol = 3, nrow = 2, width =20, height = 12) {
  if (length(x) == 0) return(NULL)
  pdf(name, width = width, height = height)
  for (i in seq(1, length(x), by = ncol*nrow)) {
    #print(names(x)[i])
    j  <- min(i + ncol*nrow - 1, length(x))
    do.call(grid.arrange, c(x[i:j], ncol = ncol, nrow = nrow))
  }
  dev.off() %>% invisible

}


## Function to generate an ordered row IDs based on the number of rows
genRowIDs <- function(x) {
  c(paste0(LETTERS,0), paste0(rep(LETTERS, each = 26), rep(LETTERS, times = 26)))[seq(1,x)]
}

## Function to generate an ordered column IDs based on the number of columns
genColIDs <- function(x) {
  formatC(seq(x),width=2,format = "d",flag="0")
}


## Function for fitting edge effect using loess local regression method
fitOneLoess <- function(plateData, useNeg, useLowConcentrations, lowConcTab, span) {

  #select "de facto" negative controll wells for surface fitting
  if (!useNeg) {
    fitIn <- plateData
  } else {
    if (! "wellType" %in% colnames(plateData)) stop("No well type information found")
    if(useLowConcentrations > 0) {
      #lowest concentrations are also used for fitting
      lowConcWell <- filter(plateData, paste0(name,"_",concentration) %in% paste0(lowConcTab$name,"_",lowConcTab$concentration))$wellID
    } else lowConcWell <- c()

    fitIn <- filter(plateData, wellType %in% "neg" | wellID %in% lowConcWell)
  }

  #fit the surface using loess model
  lo <- loess(normVal ~ numRowID + numColID, data = fitIn, control = loess.control(surface="direct"), span=span)
  plateData$edgeFactor <-  predict(lo, newdata = plateData)
  plateData
}

## Function for fitting edge effect using 2D sigmoid model (experimental feature)

fitOneSigmoid <- function(plateData, useNeg, useLowConcentrations, lowConcTab) {
  if (! "wellType" %in% colnames(plateData)) stop("No well type information found")
  #select "de facto" negative controll wells for surface fitting
  if (!useNeg) {
    fitIn <- filter(plateData, ! wellType %in% "pos")
  } else {
    if(useLowConcentrations > 0) {
      #lowest concentrations are also used for fitting
      lowConcWell <- filter(plateData, paste0(name,"_",concentration) %in% paste0(lowConcTab$name,"_",lowConcTab$concentration))$wellID
    } else lowConcWell <- c()

    fitIn <- filter(plateData, wellType %in% "neg" | wellID %in% lowConcWell)
  }


  #estimate initial values based on the values on the plate
  dr1 = median(filter(fitIn, numRowID == min(numRowID))$normVal)
  cr2 = 1- median(filter(fitIn, numRowID == max(numRowID))$normVal)
  dc1 = median(filter(fitIn, numColID == min(numColID))$normVal)
  cc2 = 1- median(filter(fitIn, numColID == max(numColID))$normVal)
  er1=max(fitIn$numRowID)*0.1
  er2=max(fitIn$numRowID) + 1 - er1
  ec1=max(fitIn$numColID)*0.1
  ec2=max(fitIn$numColID) + 1 - ec1
  sr1=0.5
  sr2=0.5
  sc1=0.5
  sc2=0.5

  paraIni <- c(dr1=dr1,cr2=cr2,dc1=dc1,cc2=cc2,er1=er1,er2=er2,ec1=ec1,ec2=ec2,sr1=sr1,sr2=sr2,sc1=sc1,sc2=sc2)

  esFun <- function(data, par) {
    with(data, sum(((normVal - (1 + (par[1]-1)/(1+exp((numRowID-par[5])/par[9])) - par[2] + par[2]/(1+exp((numRowID-par[6])/par[10])))*(1 + (par[3]-1)/(1+exp((numColID-par[7])/par[11])) - par[4] + par[4]/(1+exp((numColID-par[8])/par[12])))))^2))
  }

  predictFun <- function(numRowID,numColID, par) {
    (1 + (par[1]-1)/(1+exp((numRowID-par[5])/par[9])) - par[2] + par[2]/(1+exp((numRowID-par[6])/par[10])))*(1 + (par[3]-1)/(1+exp((numColID-par[7])/par[11])) - par[4] + par[4]/(1+exp((numColID-par[8])/par[12])))
  }

  #estimatePlate$y <- estimatePlate$normVal
  edgeFactors <- tryCatch({
    sig2D <- optim(par = paraIni, esFun, data = fitIn)
    predictFun(plateData$numRowID, plateData$numColID , sig2D$par)
  }, error = function(err) {
    warning("Edge effect estimate using 2D-sigmoid model failed. Set edgeEffect value to 1")
    rep(1, nrow(plateData))
  })

  plateData$edgeFactor <- edgeFactors
  plateData
}


#function to calculate AUC using linear-log trapezoidal method
calcAUC <- function(value, conc, viabCut=NULL) {
  valueConc <- tibble(viab=value, conc=conc) %>%
    filter(!is.na(conc) & conc > 0) %>% #remove zero concentration, not included in calculation of AUC.
    group_by(conc) %>% summarise(viab = mean(viab)) %>%
    ungroup() %>%  #make sure concentrations are unique
    arrange(conc)

  #censor the vlaues based on viability cut-off
  if (!is.null(viabCut)) {
    valueConc <- mutate(valueConc,viab = viab + 1-viabCut) %>% mutate(viab=ifelse(viab >1,1,viab))
  }

  if (nrow(valueConc) %in% c(0,1)) {
    return(mean(valueConc$viab, na.rm = TRUE))
  } else {
    areaTotal <- 0
    for (i in seq(nrow(valueConc)-1)) {
      areaEach <- (valueConc$viab[i] + valueConc$viab[i+1])*log(valueConc$conc[i+1]/valueConc$conc[i])*0.5
      areaTotal <- areaTotal + areaEach
    }
    aucRatio <- areaTotal/log(valueConc[nrow(valueConc),]$conc/valueConc[1,]$conc)
    return(aucRatio)
  }
}



# Summarise the IC50 fitting results (stand alones)
sumIC50 <- function(formula, data = NULL) {
  edgeFactors <- tryCatch({
    parm_fit <- dr4pl(formula,data)
    res <- data.frame(t(parm_fit$parameters))
  }, error = function(err) {
    warning("Curve fitting failed, NA values generated. ")
    res <- data.frame(t(structure(rep(NA,4),names=c("UpperLimit","IC50","Slope","LowerLimit"))))
  })

  return(res)
}


# Model object for fitting 4 parameter logistic (4PL) models. (Used for plotting)
fitIC50 <- function(formula, data = NULL, weights = NULL, logDose = NULL, ...) {
    if (! is.null(data) ) {
      modelFrame <- model.frame(formula, data)
    } else {
      modelFrame <- model.frame(formula)
    }

    if (!is.null(logDose)) {
      modelFrame[,2] <- logDose^modelFrame[,2]
    }

    parm_fit <- dr4pl(modelFrame[,2],modelFrame[,1])
    newModel <- list(model = modelFrame, formula = formula, parm_fit = parm_fit, logDose = logDose)

    class(newModel) <- "fitIC50"
    return(newModel)
}

#' Generic function for ic50 class generated from fitIC50 function.
predict.fitIC50 <- function(object, newdata = NULL, se.fit = FALSE, level = 0.95 ,
                            interval = c("none", "confidence", "prediction"), ...) {

  if (is.null(newdata))
    newdata <- object$model else
      newdata <- newdata

    params <- object$parm_fit$parameters
    logDose <- object$logDose

    a <- min(params[1],params[4])
    d <- max(params[4],params[1])
    c <- params[2]
    b <- params[3]

    predY <- function(x) {
      y = d + (a-d)/(1+(x/c)^b)
      names(y) <- NULL
      return(y)
    }

    if (is.vector(newdata)) {
      conc <- newdata
    } else if (is.data.frame(newdata)) {
      if (ncol(newdata) == 1) {
        conc <- newdata[,1]
      } else {
        concName <- as.character(object$formula)[3]
        conc <- newdata[,concName]
      }
    }

    if (!is.null(logDose)) conc <- logDose^conc
    res <- sapply(conc, predY)

    return(res)
}
