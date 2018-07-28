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
  sprintf("%02s",seq(x))
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
