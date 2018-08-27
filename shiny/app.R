# define modules

library(shiny)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(plotly)
library(DrugScreenExplorer)
library(DT)
load("./shinyData.RData")

ui <- navbarPage("DrugScreenExplorer", inverse = TRUE, id = "tabs",
                 tabPanel("Sample overview",
                          titlePanel("Overview of samples in the screen"),
                          dataTableOutput("sampleInfo")),
                 tabPanel("Plate plot",
                          titlePanel("Heatmap plot of drug screen plates"),
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("patientSeleBox"),
                              uiOutput("sampleSeleBox"),
                              uiOutput("plateSeleBox"),
                              uiOutput("valueSeleRadio"),
                              uiOutput("ifCensorCheck"),
                              uiOutput("censorLimitBox"),
                              uiOutput("ifCorCheck"), width =4),
                            mainPanel(tags$style(type="text/css",
                                                 ".shiny-output-error { visibility: hidden; }",
                                                 ".shiny-output-error:before { visibility: hidden; }"),
                            plotlyOutput("plateHeatmap", width = 700, height = 400),
                            column(dataTableOutput("infoTab"),width = 4),
                            column(plotOutput("linePlot", width=400, height = 300), width = 4), width = 8))),
                 tabPanel("Dose-response curves",
                          titlePanel("Plot dose-response curves for selected drugs and samples"),
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("nameDrugBox"),
                              uiOutput("seleGroup"),
                              conditionalPanel(condition = "input.selectGroup != 'All samples'"
                                               ,uiOutput("seleFeature")),
                              uiOutput("colorBox"),
                              uiOutput("ifCorCheck1"),
                              checkboxInput("ifIC50", "Fit IC50 curve", value = FALSE)
                            ),
                            mainPanel(tags$style(type="text/css",
                                                 ".shiny-output-error { visibility: hidden; }",
                                                 ".shiny-output-error:before { visibility: hidden; }"),
                              plotlyOutput("curvePlot",width = 600, height= 400)
                            ))),
                 tabPanel("Clustering",
                          titlePanel("Clustering samples or drugs based on drug response"),
                          sidebarLayout(
                            sidebarPanel(uiOutput("seleGroup1"),
                                      conditionalPanel(condition = "input.selectGroup1 != 'All samples'"
                                                       ,uiOutput("seleFeature1")),
                                      radioButtons("plotType","Plot type", c("Heatmap","PCA","t-SNE","Drug-drug correlation"),
                                                   selected = "Heatmap"),
                                      conditionalPanel(condition = "input.plotType == 't-SNE'",
                                                       numericInput("perplex", "Perplexity", value =20),
                                                       numericInput("Niter", "Iteration round", value = 1000)),
                                      conditionalPanel(condition = "input.plotType == 'Heatmap'",
                                                       radioButtons("ifRowNorm","Normalization", c("None","Row","Column"), selected = "None")),
                                      uiOutput("ifCorCheck2"),
                                      uiOutput("colorBox1"),
                                      checkboxInput("ifCensor1","Fix value range", value = FALSE),
                                      uiOutput("censorLimitBox1"),
                                      sliderInput("topPercent", "Top variant drugs", min = 10, max = 100, value = 100, step = 10),
                                      downloadButton("download", "Download as pdf"),
                                      numericInput("figWidth", label = "Figure width", value = 16),
                                      numericInput("figHeight", label = "Figure height", value = 15)
                            ),
                            mainPanel(
                              tags$style(type="text/css",
                                         ".shiny-output-error { visibility: hidden; }",
                                         ".shiny-output-error:before { visibility: hidden; }"
                              ),
                              actionButton("doPlot","Plot!", align = "center" ),
                              plotOutput("clusterPlot", width=800, height = 800)
                            ))))


server <- function(input, output, session) {

  #Decide which panel to show or hide based on input data columns
  if (! all(c("normVal","concentration","name") %in% colnames(screenData))) {
    hideTab(inputId = "tabs", target = "Dose-response curves")
  }

  ##First panel: sample over view and subsetting##

  # A table to store unique samples information
  infoTable <- screenData[,sampleAnnotations] %>%
    distinct(sampleID, .keep_all = TRUE)

  # DT table for showing the sample annotations and filtering
  output$sampleInfo <- renderDataTable({
    infoTable_new <- mutate_if(infoTable, is.character, as.factor) %>%
      select(-fileName)
    datatable(infoTable_new, filter = 'top', rownames = FALSE,
              selection = "none", style = "bootstrap")
  })

  # Reactive object for getting filtered samples
  filteredData <- reactive({
    selectedPlates <- as.character(infoTable[input$sampleInfo_rows_all,] %>% pull(fileName))
    dataNew <- filter(screenData, fileName %in% selectedPlates)
    dataNew
  })

  ##Second panel: plate heatmap for a plate##

  output$patientSeleBox <- renderUI({
    if ("patientID" %in% colnames(filteredData())) {
      selectInput("selePatient","Select a patient",
                  sort(unique(filteredData()$patientID)),
                  size=5, selectize = FALSE)
    }
  })

  output$sampleSeleBox <- renderUI({
    if (all(c("patientID","sampleID") %in% colnames(filteredData()))) {
      sampleIDs <- unique(filter(filteredData(), patientID == input$selePatient)$sampleID)
    } else if ("sampleID" %in% colnames(filteredData())) {
      sampleIDs <- unique(filteredData()$sampleID)
    } else sampleIDs <- NULL

    if (!is.null(sampleIDs)) {
      selectInput("seleSample", "Select a sample",
                  sort(sampleIDs), size=3, selectize = FALSE)
    }
  })

  output$plateSeleBox <- renderUI({
    if ("sampleID" %in% colnames(filteredData())) {
      plateIDs <- unique(filter(filteredData(), sampleID == input$seleSample)$fileName)
    } else {
      plateIDs <- unique(filteredData()$fileName)
    }

    selectInput("selePlate", "Select a plate file",
                sort(plateIDs), size=3, selectize = FALSE)
  })

  output$valueSeleRadio <- renderUI({

    if ("normVal" %in% colnames(filteredData())) {
      opts <- c( "viability", "z-score", "raw count")
    } else {
      opts <- c("z-score", "raw count")
    }
    radioButtons("ifNorm","Select value type",opts)
  })

  output$ifCorCheck <- renderUI({
    if (input$ifNorm == "viability" && ("normVal.cor" %in% colnames(filteredData()))) {
      checkboxInput("ifCorrected","Incubation effect correction", value = FALSE)
    }
  })

  output$ifCensorCheck <- renderUI({
    if(input$ifNorm %in% c("viability","z-score")) {
      checkboxInput("ifCensor","Fix value range", value = FALSE)
    }
  })

  output$censorLimitBox <- renderUI({
    if (!is.null(input$ifCensor) && input$ifCensor) {
      if (input$ifNorm == "viability"){
        upper = 1.5
        lower = 0
      } else if (input$ifNorm == "z-score") {
        upper = 4
        lower = -4
      }

      tagList(div(style="display:inline-block; width: 90px", textInput("upperLimit", "upper limit", value = upper)),
              div(style="display:inline-block; width: 90px", textInput("lowerLimit", "lower limit", value = lower)))
    }
  })

  plotTab <- reactive({
    ## Function to generate an ordered row IDs based on the number of rows
    genRowIDs <- function(x) {
      c(paste0(LETTERS,0), paste0(rep(LETTERS, each = 26), rep(LETTERS, times = 26)))[seq(1,x)]
    }
    ## Function to generate an ordered column IDs based on the number of columns
    genColIDs <- function(x) {
      sprintf("%02s",seq(x))
    }

    nCol <- length(unique(screenData$colID))
    nRow <- length(unique(screenData$rowID))
    rowSeq <- genRowIDs(nRow)
    colSeq <- genColIDs(nCol)
    matPlate <- filter(screenData, fileName == input$selePlate) %>%
      mutate(rowID = factor(rowID, levels = rev(rowSeq))) %>%
      mutate(colID = factor(colID, levels = colSeq))

    if (input$ifNorm == "viability") {
      if (!is.null(input$ifCorrected) && input$ifCorrected)
        matPlate <- mutate(matPlate, val = normVal.cor)
      else matPlate <- mutate(matPlate, val = normVal)
    } else if (input$ifNorm == "raw count") {
      matPlate <- mutate(matPlate, val = value)
    } else if (input$ifNorm == "z-score") {
      matPlate <- mutate(matPlate, val = (value-mean(value))/sd(value))
    }

    if (!is.null(input$ifCensor) && input$ifCensor) {
      matPlate <- mutate(matPlate, val = if_else(val > as.numeric(input$upperLimit),
                                           as.numeric(input$upperLimit), val)) %>%
        mutate(val = if_else(val < as.numeric(input$lowerLimit),
                                   as.numeric(input$lowerLimit), val))
    }

    matPlate
  })

  #plot the heatmap
  output$plateHeatmap <- renderPlotly({
    p <- ggplot(plotTab(), aes(x = colID, y = rowID, fill = val, label = name)) +
      geom_tile(color = "grey80") + theme_minimal() +
      theme(axis.text = element_text(size = 6), axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5),
            plot.margin = margin(5, 5, 5, 5)) +
      labs(x = "", y = "", title = input$selePlate)

    if(input$ifNorm %in% c("viability","z-score")) {
      if (!is.null(input$ifCensor) && input$ifCensor) {
        p <- p + scale_fill_gradient2(high = "red", low = "blue",
                                      mid = "white", midpoint = 1,
                                      limits = c(as.numeric(input$lowerLimit),
                                                 as.numeric(input$upperLimit)),
                                      name = input$ifNorm)
      } else {
        p <- p + scale_fill_gradient2(high = "red", low = "blue",
                                      mid = "white", midpoint = 1,
                                      name = input$ifNorm)
      }
    } else if (input$ifNorm == "raw count") {
      p <- p + scale_fill_gradient(low = "blue", high = "red",
                                   name = input$ifNorm)
    }


    ggplotly(p) %>% config(displayModeBar = F)
  })

  #a table for selected plate and well information
  output$infoTab <- renderDataTable({
    s <- event_data("plotly_click")

    #Add potential interesting columns that may or may not exist
    addCol <- intersect(c("name","wellType","concentration"), colnames(plotTab()))
    showAnno <- unique(c(sampleAnnotations, addCol))
    showAnno <- showAnno[showAnno != "fileName"]
    infoTab <- filter(plotTab(), rowID == levels(plotTab()$rowID)[s$y],
                      colID == levels(plotTab()$colID)[s$x]) %>%
      select(!!showAnno) %>%
      mutate(`plate median value` = median(plotTab()$value)) %>%
      t

    datatable(infoTab, colnames = "",selection = 'none', options = list(dom = "t"))})

  #prepare table for the dose-response plot on the first panel
  linePlotTab1 <- reactive({
    if (any(!c("name","concentration") %in% colnames(plotTab()))) {
      NULL
    } else {
      s <- event_data("plotly_click")
      drugName <- unique(filter(plotTab(), rowID == levels(plotTab()$rowID)[s$y],
                        colID == levels(plotTab()$colID)[s$x])$name)
      plotTab <- filter(plotTab(), name == drugName)
      plotTab
    }
  })

  #dose-response curve for selected drugs if possible.
  output$linePlot <- renderPlot({

    if (any(!c("name","concentration") %in% colnames(plotTab()))) {
      text <- "drug name, concentration and normalized viability\n must be available in order to show\n dose-response curve"
      p <- ggplot() +
        annotate("text", x = 4, y = 25, size=3, label = text, col = "red") +
        theme_void() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
    } else {
      s <- event_data("plotly_click")
      drugName <- unique(filter(plotTab(), rowID == levels(plotTab()$rowID)[s$y],
                                colID == levels(plotTab()$colID)[s$x])$name)
      plotTab <- filter(plotTab(), name == drugName)
      p <- ggplot(plotTab, aes(x = concentration, y = val)) +
        geom_point(col = "red") + geom_line(col = "red") + scale_x_log10() +
        theme_bw() + ylab(input$ifNorm) + ggtitle(drugName) +
        theme(plot.title = element_text(hjust=0.5, face = "bold", color = "red"),
              axis.text = element_text(size = 13),
              axis.title = element_text(size=13))
      if (input$ifCensor) {
        p <- p + coord_cartesian(ylim = c(as.numeric(input$lowerLimit),
                                          as.numeric(input$upperLimit)))
      }
    }
    p
  })

  ## Third panel: dose-response curves ##

  #select a drug name
  output$nameDrugBox <- renderUI({
    allDrugs <- filter(filteredData(), wellType == "sample") %>%
      distinct(name) %>% pull(name)
    selectInput("nameDrug","Select a drug",sort(allDrugs))
  })

  #select patient group for subsetting
  output$seleGroup <- renderUI({
    selectInput("selectGroup", "Subset samples by", c("All samples", sampleAnnotations))
  })

  #select features for grouping of the patients
  output$seleFeature <- renderUI({
    if (input$selectGroup != "All samples") {
      allGroups <- unique(pull(filteredData(),input$selectGroup))
      selectInput("seleSub","Groups included",allGroups,
                  size = 5, selectize = FALSE, multiple = TRUE,
                  selected = allGroups[1])
    }
  })

  #select features for color scheme
  output$colorBox <- renderUI({
    selectInput("seleColor","Color by",sampleAnnotations,
                size = 5, selectize = FALSE, multiple = FALSE,
                selected = sampleAnnotations[1])
  })

  #show option for incubation effect correction if "normVal.cor" is present
  output$ifCorCheck1 <- renderUI({
    if ("normVal.cor" %in% colnames(filteredData())) {
      checkboxInput("ifCorrected1","Incubation effect correction", value = FALSE)
    }
  })



  #data table for dose-response plot
  doseTab <- reactive({

    #subset the screen data according to the selection
    if (input$selectGroup == "All samples") {
      plotTab <- filter(filteredData(), name == input$nameDrug)
    } else {
      selePlates <- filteredData()[unlist(filteredData()[,input$selectGroup]) %in% input$seleSub,]$fileName
      plotTab <- filter(filteredData(), fileName %in% selePlates, name == input$nameDrug)
    }

    plotTab
  })

  #reactive object to make the dose-response curves
  doseCurves <- reactive({

      fmt_dcimals <- function(decimals=0){
        # return a function responpsible for formatting the
        # axis labels with a given number of decimals
        function(x) as.character(round(x,decimals))
      }

      plotTab <- doseTab()

      if (!is.null(input$ifCorrected1) && input$ifCorrected1) {
        valueType <- "normVal.cor"
      } else {
        valueType <- "normVal"
      }

      p <- ggplot(data=plotTab, aes_string(x="concentration", y=valueType,
                                              col = input$seleColor, group = "sampleID")) +
        geom_point(pch=16, size=4) + theme_classic() + scale_x_log10(labels = fmt_dcimals(3)) +
        theme(axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              plot.title = element_text(hjust = 0.5, face = "bold", color = "red")) +
        ylab("Viability") + xlab("Concentration") + coord_cartesian(ylim = c(0,1.5)) + ggtitle(input$nameDrug)
      if (input$ifIC50) p <- p + geom_smooth(method = "fitIC50", se=FALSE, method.args= list(logDose =10)) else
        p <- p + stat_summary(aes(group=sampleID), fun.y=mean, geom="line")

      p

  })

  output$curvePlot <- renderPlotly({
    ggplotly(doseCurves())  %>% config(displayModeBar = F)
  })

  ## Fourth panel:  Clustering plot ##

  #select patient group for subsetting
  output$seleGroup1 <- renderUI({
    selectInput("selectGroup1", "Subset samples by", c("All samples", sampleAnnotations))
  })

  #select features for grouping of the patients
  output$seleFeature1 <- renderUI({
    if (input$selectGroup1 != "All samples") {
      allGroups <- unique(pull(filteredData(),input$selectGroup1))
      selectInput("seleSub1","Groups included",allGroups,
                  size = 5, selectize = FALSE, multiple = TRUE,
                  selected = allGroups[1])
    }
  })


  #select features for color scheme on the fourth pannel (for clustering plot)
  output$colorBox1 <- renderUI({
      selectInput("seleColor1","Color by",sampleAnnotations,
                  size = 5, selectize = FALSE, multiple = FALSE,
                  selected = sampleAnnotations[1])
  })

  #show option for incubation effect correction if "normVal.cor" is present
  output$ifCorCheck2 <- renderUI({
    if ("normVal.cor" %in% colnames(filteredData())) {
      checkboxInput("ifCorrected2","Incubation effect correction", value = FALSE)
    }
  })


  output$censorLimitBox1 <- renderUI({
    if (!is.null(input$ifCensor1) && input$ifCensor1) {
    upper = 1.5
    lower = 0
    tagList(div(style="display:inline-block; width: 90px", textInput("upperLimit1", "upper limit", value = upper)),
            div(style="display:inline-block; width: 90px", textInput("lowerLimit1", "lower limit", value = lower)))
    }
  })


  #table for the clustering plot on panel 3
  plotMat <- reactive({

    #subset the screen data according to the selection
    if (input$selectGroup1 == "All samples") {
      plotTab <- filter(filteredData(), wellType == "sample")
    } else {
      selePlates <- filteredData()[unlist(filteredData()[,input$selectGroup1]) %in% input$seleSub1,]$fileName
      plotTab <- filter(filteredData(), fileName %in% selePlates, wellType == "sample")
    }

    if (!is.na(input$ifCorrected2) & input$ifCorrected2) {
      plotTab <- mutate(plotTab, viab = normVal.cor_auc)
    } else {
      plotTab <- mutate(plotTab, viab = normVal_auc)
    }

    if(input$ifCensor1) {
      plotTab <- mutate(viab = ifelse(viab > as.numeric(input$upperLimit1), as.numeric(input$upperLimit1),viab)) %>%
        mutate(viab = ifelse(viab < as.numeric(input$lowerLimit1), as.numeric(input$lowerLimit1),viab))
    }

    viabMat <- group_by(plotTab, name, sampleID) %>%
      summarise(viab = mean(viab))
    idMap <- tibble(origin = as.character(viabMat$sampleID))
    viabMat <- spread(viabMat, key = sampleID, value = viab) %>%
      data.frame() %>% column_to_rownames("name") %>%
      as.matrix()
    #in case the column names begins with number and changed by R
    idMap$new <- colnames(viabMat)

    viabMat <- viabMat[complete.cases(viabMat),] #only use drugs that present in all samples

    #variant filtering
    sds <- apply(viabMat,1,sd)
    viabMat <- viabMat[sds >= quantile(sds, 1 - input$topPercent/100),]


    list(viabMat = viabMat, idMap = idMap)
  })

  #reactive object for PCA
  tabPCA <- reactive({
    viabMat <- plotMat()$viabMat
    pcaRes <- prcomp(t(viabMat), scale. = FALSE, center = TRUE)
    varExp <- (pcaRes$sdev^2 / sum(pcaRes$sdev^2))*100
    pcaTab <- data.frame(pcaRes$x[,c(1,2)])
    names(varExp) <- colnames(pcaRes$x)
    rownames(pcaTab) <- colnames(viabMat)
    list(pcaTab = pcaTab,varExp = varExp)
  })

  #reactive object for tSNE
  tabTSNE <- reactive({
    library(Rtsne)
    viabMat <- plotMat()$viabMat
    #prepare distance matrix
    distViab <- dist(t(viabMat))
    tsneRes <- Rtsne(distViab, perplexity = input$perplex, theta = 0,
                     max_iter = input$Niter, is_distance = TRUE, dims =2)
    tsneRes <- tsneRes$Y
    rownames(tsneRes) <- labels(distMat)
    colnames(tsneRes) <- c("x","y")
    plotTab <- data.frame(tsneRes)
    rownames(plotTab) <- colnames(viabMat)
    plotTab
  })

  #prepare the clustering plot on panel 3
  plotCluster <- eventReactive(input$doPlot, {
    viabData <- filteredData()
    viabMat <- plotMat()$viabMat
    idMap <- plotMat()$idMap
    annoCol <- tibble(sampleID.new = colnames(viabMat),
      sampleID = idMap[match(colnames(viabMat), idMap$new),]$origin) %>%
      left_join(select(viabData, sampleID, input$seleColor1), by = "sampleID") %>%
      select(-sampleID) %>% distinct(sampleID.new, .keep_all = TRUE) %>%
      data.frame() %>% column_to_rownames("sampleID.new")

    if (input$plotType == "Heatmap"){

      #perform hierachical clustering on columns

      hc <- hclust(as.dist(1-cor(viabMat)), method = "ward.D2")
      hr <- hclust(as.dist(1-cor(t(viabMat))), method = "ward.D2")
      #hc <- hclust(dist(t(viabMat), method = "euclidean"), method = "ward.D2")

      #perfrom hierachical clustering on rows, normalized or not
      if (input$ifRowNorm == "Row") {
        viabMat <- t(scale(t(viabMat)))
      } else if (input$ifRowNorm == "Column") {
        viabMat <- scale(viabMat)
      }

      annoCol <- annoCol[, colSums(!is.na(annoCol)) != 0, drop = FALSE]

      if(ncol(annoCol) != 0) {
        g <- pheatmap(viabMat,color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100), scale = "none",
                      annotation_col = annoCol, cluster_rows = hr,
                      cluster_cols = hc)$gtable
      } else {
        #if nothing to annotate
        g <- pheatmap(viabMat,color = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100), scale = "none",
                      cluster_rows = hr, cluster_cols = hc)$gtable
      }
      g

    } else if (input$plotType == "PCA") {
      pcaTab <- cbind(tabPCA()$pcaTab,annoCol)
      varExp <- tabPCA()$varExp
      pcaTab$sampleID <- rownames(pcaTab)
      g <- ggplot(pcaTab, aes_string(x="PC1",y="PC2",color = input$seleColor1)) + geom_point() + theme_bw() +
        xlab(sprintf("PC1 (%2.1f%%)",varExp[1])) + ylab(sprintf("PC2 (%2.1f%%)",varExp[2]))
      g

    } else if (input$plotType == "t-SNE") {
      plotTab <- cbind(tabTSNE(),annoCol)
      plotTab$sampleID <- rownames(plotTab)
      g <- ggplot(plotTab, aes_string(x="x",y="y",color = input$seleColor1)) + geom_point() + theme_bw()
      g

    } else if (input$plotType == "Drug-drug correlation") {
      g <- pheatmap(cor(t(viabMat),method="spearman"))
      g$gtable
    }
  })

  #clustering plot on panel 3
  output$clusterPlot <- renderPlot({
    plotCluster()
  })

  #download the cluster plot
  output$download <- downloadHandler(
    filename = paste0(input$plotType, '.pdf', sep=''),
    content = function(file) {
      ggsave(file, plot = plotCluster(),
             device = "pdf", width = input$figWidth, height = input$figHeight,
             limitsize = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
