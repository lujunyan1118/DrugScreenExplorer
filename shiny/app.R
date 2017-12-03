# define modules

library(tidyverse)
library(DT)
load("./shinyData.RData")

subsetSampleUI <- function(id, screenData, sampleAnnotations) {
  ns <- NS(id)
  tagList(radioButtons(ns("ifSubset"),"Sample selection",
                       c("all samples","subset samples"),
                       selected = "all samples"),
          conditionalPanel(paste0("input['", ns("ifSubset"), "'] == 'subset samples'"),
                           selectInput(ns("subsetBy"), "Subset samples by",
                                       choices = sampleAnnotations),
                           uiOutput(ns("subsetSele"))))
}

seleSample <- function(input, output, session, screenData){
   input$subsetSele <- renderUI({
     factors <- pull(screenData, input$subsetBy)
     if (is.character(factors) | is.factor(factors)) {
       #the annotation type is
     }
   })

}



ui <- fluidPage(navbarPage("DrugScreenExplorer", inverse = TRUE,
                           tabPanel("Sample overview",
                                    titlePanel("Overview of samples in the screen")
                                    )))


server <- function(input, output, session) {}

shinyApp(ui = ui, server = server)
