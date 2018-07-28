# define modules

library(shiny)
library(tidyverse)
library(plotly)
library(DrugScreenExplorer)
library(DT)
load("../vignettes/src/shiny/shinyData.RData")

ui <- navbarPage("DrugScreenExplorer", inverse = TRUE,
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
                              uiOutput("ifCorCheck1"), width =4),
                            mainPanel(tags$style(type="text/css",
                                                 ".shiny-output-error { visibility: hidden; }",
                                                 ".shiny-output-error:before { visibility: hidden; }"),
                            plotlyOutput("plateHeatmap", width = 700, height = 400),
                            column(dataTableOutput("infoTab"),width = 4),
                            column(plotOutput("linePlot", width=400, height = 300), width = 4), width = 8))))


server <- function(input, output, session) {

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

  #a test output
  #output$testOut <- renderTable(filteredData())

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

  output$ifCorCheck1 <- renderUI({
    if (input$ifNorm == "viability" && ("normVal.cor" %in% colnames(filteredData()))) {
      checkboxInput("ifCorrected1","Incubation effect correction", value = FALSE)
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
      if (!is.null(input$ifCorrected1) && input$ifCorrected1)
        matPlate <- mutate(matPlate, val = normVal.cor)
      else matPlate <- mutate(matPlate, val = normVal)
    } else if (input$ifNorm == "raw count") {
      matPlate <- mutate(matPlate, val = value)
    } else if (input$ifNorm == "z-score") {
      matPlate <- mutate(matPlate, val = (value-mean(value))/sd(value))
    }

    if (!is.null(input$ifCensor) & input$ifCensor) {
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
      if (!is.null(input$ifCensor) & input$ifCensor) {
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

}

shinyApp(ui = ui, server = server)
