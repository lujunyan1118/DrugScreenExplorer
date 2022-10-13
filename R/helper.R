# Some helper functions

## Function to read screen data from one plate
readPlate <- function(plateFile, rowRange = c(3, 18), colRange = 2, sep = "[;\t]", commaAsDecimal = FALSE,
                      allowNA = FALSE) {

  if (grepl(",", sep) & commaAsDecimal)
    stop("Comma can not be used as a separator and the demical at the same time.")

  txt <- readLines(plateFile)

  if (length(rowRange) == 1) {
    rowRange <- c(rowRange, length(txt))
  }

  sp <- lapply(seq(rowRange[1], rowRange[2]), function(x) {
    if (commaAsDecimal) {
      txtSp <- gsub("[,]", ".", txt[x])
    } else {
      txtSp <- txt[x]  #only ',' or '.' as decimal
    }
    txtSp <- strsplit(txtSp, split = sep)[[1]]
    #txtSp <- txtSp[txtSp != ""]
    if (length(colRange) == 1)
      colRange <- c(colRange, length(txtSp))
    txtSp <- txtSp[seq(colRange[1], colRange[2])]
  })

  # check whether all the rows have the same length.
  rowNum <- length(sp)
  colNum <- unique(vapply(sp, length, numeric(1)))

  if (length(colNum) != 1) {
    stop("Numbers of the values in rows are different.")
  }

  rowID <- genRowIDs(rowNum)
  colID <- genColIDs(colNum)
  rowID <- rep(rowID, each = colNum)
  colID <- rep(colID, times = rowNum)

  df <- tibble::tibble(wellID = paste0(rowID, colID), rowID = rowID,
                       colID = colID, value = as.numeric(unlist(sp)))

  if (any(is.na(df$value))) {
    if (allowNA) {
      warning("NA values detected in plate")
    } else {
      stop("NA values detected in plate. If NA values should be allowed, please set the 'allowNA = TRUE'")
    }
  }

  return(df)
}



## Function to make pdf for a list of ggplot objects
makepdf <- function(x, name, ncol = 3, nrow = 2, width = 20, height = 12) {

    if (length(x) == 0)
        return(NULL)
    pdf(name, width = width, height = height)
    for (i in seq(1, length(x), by = ncol * nrow)) {
        # print(names(x)[i])
        j <- min(i + ncol * nrow - 1, length(x))
        do.call(gridExtra::grid.arrange, c(x[seq(i, j)], ncol = ncol, nrow = nrow))
    }

    dev.off() %>%
        invisible

}


## Function to generate an ordered row IDs based on the number of rows
genRowIDs <- function(x) {
    if (x > 702)
        stop("Rows exceeds maximal supported number (n=702)")

    allID <- c(paste0(LETTERS, 0), paste0(rep(LETTERS, each = 26), rep(LETTERS, times = 26)))

    return(allID[seq(1, x)])
}

## Function to generate an ordered column IDs based on the number of columns
genColIDs <- function(x) {
    formatC(seq(x), width = 2, format = "d", flag = "0")
}


## Function for fitting edge effect using loess local regression method
fitOneLoess <- function(plateData, useNeg, useLowConcentrations, lowConcTab, span) {

    # select 'de facto' negative controll wells for surface fitting
    if (!useNeg) {
        fitIn <- plateData
    } else {

        if (!"wellType" %in% colnames(plateData))
            stop("No well type information found")
        if (useLowConcentrations > 0) {

            # lowest concentrations are also used for fitting
            lowConc <- paste0(lowConcTab$name, "_", lowConcTab$concentration)
            lowConcWell <- plateData[paste0(plateData$name, "_", plateData$concentration) %in% lowConc, ]$wellID

        } else lowConcWell <- c()

        fitIn <- plateData[plateData$wellType %in% "neg" | plateData$wellID %in% lowConcWell, ]
    }

    # fit the surface using loess model
    lo <- loess(normVal ~ numRowID + numColID, data = fitIn, control = loess.control(surface = "direct"), span = span)
    plateData$edgeFactor <- predict(lo, newdata = plateData)
    plateData
}

# Get initial values based on the values on the plate
getIniParm <- function(fitIn) {
    # estimate initial values based on the values on the plate
    dr1 <- median(fitIn[fitIn$numRowID == min(fitIn$numRowID), ]$normVal)
    cr2 <- 1 - median(fitIn[fitIn$numRowID == max(fitIn$numRowID), ]$normVal)
    dc1 <- median(fitIn[fitIn$numColID == min(fitIn$numColID), ]$normVal)
    cc2 <- 1 - median(fitIn[fitIn$numColID == max(fitIn$numColID), ]$normVal)

    er1 <- max(fitIn$numRowID) * 0.1
    er2 <- max(fitIn$numRowID) + 1 - er1
    ec1 <- max(fitIn$numColID) * 0.1
    ec2 <- max(fitIn$numColID) + 1 - ec1
    sr1 <- 0.5
    sr2 <- 0.5
    sc1 <- 0.5
    sc2 <- 0.5

    paraIni <- c(dr1 = dr1, cr2 = cr2, dc1 = dc1,
                 cc2 = cc2, er1 = er1, er2 = er2,
                 ec1 = ec1, ec2 = ec2, sr1 = sr1, sr2 = sr2,
                 sc1 = sc1, sc2 = sc2)

    return(paraIni)
}

#function for 2D sigmoid curv-fitting
esFun <- function(data, par) {
    with(data, sum(((normVal - (1 + (par[1] - 1)/(1 + exp((numRowID - par[5])/par[9])) -
                                    par[2] + par[2]/(1 + exp((numRowID - par[6])/par[10]))) *
                         (1 + (par[3] - 1)/(1 + exp((numColID - par[7])/par[11])) -
                              par[4] + par[4]/(1 + exp((numColID - par[8])/par[12])))))^2))
}

#function for predict values based on 2D sigmoid model
predictFun <- function(numRowID, numColID, par) {
    (1 + (par[1] - 1)/(1 + exp((numRowID - par[5])/par[9])) - par[2] + par[2]/(1 + exp((numRowID - par[6])/par[10]))) *
        (1 + (par[3] - 1)/(1 + exp((numColID - par[7])/par[11])) - par[4] + par[4]/(1 + exp((numColID - par[8])/par[12])))
}

## Function for fitting edge effect using 2D sigmoid model (experimental feature)
fitOneSigmoid <- function(plateData, useNeg, useLowConcentrations, lowConcTab) {

    if (!"wellType" %in% colnames(plateData))
        stop("No well type information found")

    # select 'de facto' negative control wells for surface fitting
    if (!useNeg) {
        fitIn <- plateData[!plateData$wellType %in% "pos", ]  #do not use positive controls

    } else {
        if (useLowConcentrations > 0) {
            # lowest concentrations are also used for fitting
            lowConc <- paste0(lowConcTab$name, "_", lowConcTab$concentration)
            lowConcWell <- plateData[paste0(plateData$name, "_", plateData$concentration) %in% lowConc, ]$wellID

        } else lowConcWell <- c()

        fitIn <- plateData[plateData$wellType %in% "neg" | plateData$wellID %in% lowConcWell, ]
    }

    paraIni <- getIniParm(fitIn)

    edgeFactors <- tryCatch({
        sig2D <- optim(par = paraIni, esFun, data = fitIn)
        predictFun(plateData$numRowID, plateData$numColID, sig2D$par)
    }, error = function(err) {
        warning("Edge effect estimate using 2D-sigmoid model failed. Set edgeEffect value to 1")
        rep(1, nrow(plateData))
    })

    plateData$edgeFactor <- edgeFactors
    plateData
}


# function to calculate AUC using linear-log trapezoidal method
calcAUC <- function(value, conc, viabCut = NULL) {

    valueConc <- data.frame(viab = value, conc = conc) %>%
        dplyr::filter(!is.na(conc) & conc > 0) %>%
        dplyr::group_by(conc) %>%
        dplyr::summarise(viab = mean(viab)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(conc)

    # censor the values based on viability cut-off
    if (!is.null(viabCut)) {
        valueConc <- dplyr::mutate(valueConc, viab = viab + 1 - viabCut) %>%
            dplyr::mutate(viab = ifelse(viab > 1, 1, viab))
    }

    if (nrow(valueConc) %in% c(0, 1)) {
        return(mean(valueConc$viab, na.rm = TRUE))

    } else {
        areaTotal <- 0
        for (i in seq(nrow(valueConc) - 1)) {
            areaEach <- (valueConc$viab[i] + valueConc$viab[i + 1]) * log(valueConc$conc[i + 1]/valueConc$conc[i]) * 0.5
            areaTotal <- areaTotal + areaEach
        }
        aucRatio <- areaTotal/log(valueConc[nrow(valueConc), ]$conc/valueConc[1, ]$conc)
        return(aucRatio)
    }
}

# Function to correct edge effect using a specialized linear model, to fix issue of over correcting of bliss model. This
# is only an experimental funciton.
corEdgeLinear <- function(x, m) {

    inTab <- data.frame(normVal = x, edgeFactor = m)

    s <- apply(inTab, 1, function(eachRow) {
        if (eachRow[1] < 0) {
            eachRow[1]
        } else if (eachRow[1] >= 0 & eachRow[1] <= eachRow[2]) {
            eachRow[1]/eachRow[2]
        } else {
            eachRow[1] + 1 - eachRow[2]
        }
    })
    s
}

# estimate the edge effect on one plate

estimateOnePlate <- function(plateData, nLayer, onlyNeg) {
    if (nLayer < 1) {
        stop("At least only layer needs to be specified as edge")
    } else {
        nCol <- length(unique(plateData$colID))
        nRow <- length(unique(plateData$rowID))
        edgeCol <- genColIDs(nCol)[c(seq(1, nLayer), seq(nCol - nLayer + 1, nCol))]
        edgeRow <- genRowIDs(nRow)[c(seq(1, nLayer), seq(nRow - nLayer + 1, nRow))]
    }

    edgeWell <- dplyr::filter(plateData, colID %in% edgeCol, rowID %in% edgeRow) %>%
        dplyr::select(wellID, value, wellType)
    innerWell <- dplyr::filter(plateData, !wellID %in% edgeWell$wellID)

    if (onlyNeg) {
        # only use negative control wells, check wether they are present in edge and innel wells
        if (!"neg" %in% edgeWell$wellType)
            stop("No negative control wells found on the edge.")
        if (!"neg" %in% innerWell$wellType)
            stop("No negative control wells found on the inner plate")
        edgeFac <- median(edgeWell$value, na.rm=TRUE)/median(innerWell$value, na.rm = TRUE)
    } else {
        edgeFac <- median(edgeWell$value, na.rm=TRUE)/median(innerWell$value, na.rm = TRUE)
    }

    edgeFac
}

# Determin the type of input vector
detectClass <- function(x) {
    # change potential NA records to NA
    x[x %in% c("", "NA")] <- NA
    x <- na.omit(x)
    nx <- length(unique(x))

    if (nx <= 1) {
        return("nd")
    } else if (nx == 2) {
        return("binary")
    } else {
        x.num <- as.numeric(as.character(x))
        if (all(!is.na(x.num))) {
            return("continuous")
        } else {
            if (nx < length(x)) {
                return("categorical")
            } else return("nd")
        }
    }
}


#plot viability heatmap for each plates
plotViability <- function(matPlate, plateList, limits, ifCorrected) {
    # plot normalized ATP count (viability))
    if (!"normVal" %in% colnames(matPlate))
        stop("No relative viability values found, please normalize the plate first")

    if (is.null(limits)) {
        limits <- c(0, 2)
    }

    matPlate <- dplyr::mutate(matPlate, normVal = ifelse(normVal < limits[1], limits[1],
                                                         ifelse(normVal > limits[2], limits[2], normVal)))

    if (ifCorrected) {
        if (!"normVal.cor" %in% colnames(matPlate)) {
            stop("No edge-corrected viability found, please run edge effect correction first")
        } else {
            matPlate <- dplyr::mutate(matPlate, normVal = normVal.cor)
        }
    }

    plotList <- lapply(plateList, function(plateName) {
        p <- dplyr::filter(matPlate, fileName == plateName) %>%
            ggplot(aes(x = colID, y = rowID, fill = normVal)) + geom_tile(color = "grey80") +
            labs(x = "", y = "", fill = "viability", title = plateName) + theme_void() +
            theme(axis.text = element_text(size = 6), axis.ticks = element_blank(),
                  plot.title = element_text(hjust = 0.5), plot.margin = margin(5, 5,
                                                                               5, 5))
        if (!is.null(limits)) {
            p <- p + scale_fill_gradient2(limits = limits, high = "red", low = "blue",
                                          mid = "white", midpoint = 1)
        } else {
            p <- p + scale_fill_gradient2(high = "red", low = "blue", mid = "white",
                                          midpoint = 1)
        }
    })
    return(plotList)
}

#plot zscore heatmap for each plate
plotZscore <- function(matPlate, plateList, limits) {
    # calculate per-plate z-score
    calcZ <- function(x) dplyr::mutate(x, zscore = (value - mean(value))/sd(value))
    matPlate <- group_by(matPlate, fileName) %>%
        do(calcZ(.)) %>%
        ungroup()

    if (is.null(limits)) {
        limits <- c(-3, 3)
    }

    # censoring value based on limit
    matPlate <- dplyr::mutate(matPlate, zscore = ifelse(zscore < limits[1], limits[1],
                                                        ifelse(zscore > limits[2], limits[2], zscore)))

    plotList <- lapply(plateList, function(plateName) {
        p <- dplyr::filter(matPlate, fileName == plateName) %>%
            ggplot(aes(x = colID, y = rowID, fill = zscore)) + geom_tile(color = "grey80") +
            labs(x = "", y = "", fill = "Z score", title = plateName) + theme_void() +
            theme(axis.text = element_text(size = 6), axis.ticks = element_blank(),
                  plot.title = element_text(hjust = 0.5), plot.margin = margin(5, 5,
                                                                               5, 5))
        if (!is.null(limits)) {
            p <- p + scale_fill_gradient2(limits = limits, high = "red", low = "blue",
                                          mid = "white", midpoint = 0)
        } else {
            p <- p + scale_fill_gradient2(high = "red", low = "blue", mid = "white",
                                          midpoint = 0)
        }
    })

    return(plotList)
}

#plot edge effect for each plate
plotEdgeEffect <- function(matPlate, plateList, limits) {
    # only plot estiamted edge effect
    if (!"edgeFactor" %in% colnames(matPlate))
        stop("No edge effect information found, please fit edge effect first")

    if (is.null(limits)) {
        limits <- c(0, 2)
    }

    matPlate <- dplyr::mutate(matPlate, normVal = ifelse(normVal < limits[1], limits[1],
                                                         ifelse(normVal > limits[2], limits[2], normVal)))

    plotList <- lapply(plateList, function(plateName) {
        p <- dplyr::filter(matPlate, fileName == plateName) %>%
            ggplot(aes(x = colID, y = rowID, fill = edgeFactor)) + geom_tile(color = "grey80") +
            labs(x = "", y = "", fill = "viability", title = plateName) + theme_void() +
            theme(axis.text = element_text(size = 6), axis.ticks = element_blank(),
                  plot.title = element_text(hjust = 0.5), plot.margin = margin(5, 5,
                                                                               5, 5))
        if (!is.null(limits)) {
            p <- p + scale_fill_gradient2(limits = limits, high = "red", low = "blue",
                                          mid = "white", midpoint = 1)
        } else {
            p <- p + scale_fill_gradient2(high = "red", low = "blue", mid = "white",
                                          midpoint = 1)
        }
    })

    return(plotList)
}


