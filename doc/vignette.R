## ---- message=FALSE, warning=FALSE--------------------------------------------
library(DrugScreenExplorer)
library(tidyr)
library(dplyr)
library(gridExtra)
library(tools)

## -----------------------------------------------------------------------------
system.file("testData/rawData", package = "DrugScreenExplorer")

## -----------------------------------------------------------------------------
system.file("testData/plateAnno/plateAnno.csv", package = "DrugScreenExplorer")

## -----------------------------------------------------------------------------
system.file("testData/wellAnno/wellAnno.csv", package = "DrugScreenExplorer")

## -----------------------------------------------------------------------------
rawFolder <- system.file("testData/rawData", package = "DrugScreenExplorer")

## ---- eval=FALSE--------------------------------------------------------------
#  createPlateInput(rawDir = rawFolder,
#                   file = "plateAnno.csv",
#                   entries = c("sampleID", "patientID"))

## ---- eval=FALSE--------------------------------------------------------------
#  rawFolder <- system.file("testData/rawData", package = "DrugScreenExplorer")
#  createPlateInput(rawDir = rawFolder,
#                   file = "plateAnno.csv",
#                   entries = c("sampleID", "patientID"),
#                   batchAsFolder = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  createWellInput(file = "wellAnno.csv", colNum = 24, rowNum = 16,
#                  entries = c("name", "concentration"))

## ---- eval=FALSE--------------------------------------------------------------
#  createWellInput(file = "wellAnno.csv", colNum = 24, rowNum = 16,
#                  entries = c("name", "concentration"),
#                  platePerSample = 2)

## ---- message=FALSE-----------------------------------------------------------
plateFile <- system.file("testData/plateAnno/plateAnno.csv",
                         package = "DrugScreenExplorer")
wellFile <- system.file("testData/wellAnno/wellAnno.csv",
                         package = "DrugScreenExplorer")
screenData <- readScreen(rawDir = rawFolder, plateAnnotationFile = plateFile,
                         wellAnnotationFile = wellFile, 
                         rowRange = c(3,18), colRange = 2, 
                         sep = "[;\t]",
                         negWell <- c("DMSO","PBS"), posWell = c())
head(screenData)

## ----message=FALSE------------------------------------------------------------
screenData <- readScreen(rawDir = rawFolder, plateAnnotationFile = plateFile,
                         wellAnnotationFile = wellFile, 
                         rowRange = c(3,18), colRange = 2, 
                         sep = "[;\t]",
                         negWell <- c("DMSO","PBS"), posWell = c(),
                         normalization = TRUE, method = "negatives", discardLayer = 2)
head(screenData %>% dplyr::select(wellID, value, sampleID, name, concentration, wellType, normVal))

## -----------------------------------------------------------------------------
#disable normalization by setting "normalization = FALSE"
screenData <- readScreen(rawDir = rawFolder, plateAnnotationFile = plateFile,
                         wellAnnotationFile = wellFile, 
                         rowRange = c(3,18), colRange = 2, 
                         sep = "[;\t]",
                         negWell <- c("DMSO","PBS"), posWell = c(),
                         normalization = FALSE)

## -----------------------------------------------------------------------------
screenData <- normalizePlate(screenData, method = "negatives", discardLayer = 2)

head(screenData %>% dplyr::select(wellID, value, sampleID, name, concentration, wellType, normVal))

## ---- fig.height=10, fig.width=6----------------------------------------------
g <- plotRawCount(screenData, ifLog10=TRUE)
plot(g)

## ---- fig.height=24, fig.width=12---------------------------------------------
g <- plotTypeDist(screenData, ifLog10 = TRUE)
grid.arrange(grobs = g, ncol=3)

## ---- fig.width=12, fig.height=15---------------------------------------------
g <- platePlot(screenData = screenData, plotPlate = "all", plotType = "viability")
grid.arrange(grobs = g[1:8], ncol = 2)

## ---- fig.width=12, fig.height=15---------------------------------------------
g <- platePlot(screenData = screenData, plotPlate = "all", plotType = "zscore")
grid.arrange(grobs = g[1:8], ncol = 2)

## ---- fig.height=4, fig.width=6-----------------------------------------------
g <- platePlot(screenData = screenData, plotPlate = "all", plotType = "layout")
g[[1]]

## ---- eval = FALSE------------------------------------------------------------
#  makeReport(screenData = screenData, showCode = FALSE,
#             title = "Report for my drug screening project",
#             author = "Junyan Lu", ifPlatePlot = TRUE)

## ---- fig.width=10, fig.height=3----------------------------------------------
plotPlates <- c("CTGLuminescence_P0024_14PB0550",
                "CTGLuminescence_P0014_15PB0031")
subScreen <- dplyr::filter(screenData, fileName %in% plotPlates)
g <- platePlot(subScreen, plotPlate = "all", plotType = "viability")
grid.arrange(grobs=g, ncol=2)

## -----------------------------------------------------------------------------
subScreen <- fitEdgeEffect(subScreen, method = "loess", useNeg = TRUE,
                           useLowConcentrations = 1, span = 1)

## ---- fig.width=10, fig.height=3----------------------------------------------
g <- platePlot(subScreen, plotType = "edgeEffect")
grid.arrange(grobs = g, ncol =2)

## -----------------------------------------------------------------------------
subScreen <- correctEdgeEffect(subScreen, correctMethod = "bliss")

## ---- fig.width=10, fig.height=3----------------------------------------------
g <- platePlot(subScreen, plotType = "viability", ifCorrected = TRUE)
grid.arrange(grobs = g, ncol =2)

## -----------------------------------------------------------------------------
screenData <- correctEdgeEffect(screenData, correctMethod = "bliss")

## -----------------------------------------------------------------------------
screenData <- summariseScreen(screenData, method = c("average","AUC"))

## ---- eval=FALSE--------------------------------------------------------------
#  makeShiny(screenData, sampleAnnotations = c("batch","sampleID","patientID","Mutation","Expression","Group"))

