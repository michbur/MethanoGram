
library(shiny)
library(DT)
library(biogram)
library(AmyloGram)
library(mlr)
library(dplyr)

options(shiny.maxRequestSize=10*1024^2)

options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 15
))

my_DT <- function(x)
  datatable(x, escape = FALSE, extensions = 'Buttons',
            filter = "none", rownames = FALSE)

source("functions.R")
load("pred_list.RData")

shinyServer(function(input, output) {
  
  prediction <- reactive({
    
    if (!is.null(input[["seq_file"]]))
      input_sequences <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({
      if (!is.null(input[["text_area"]]))
        if(input[["text_area"]] != "")
          input_sequences <- read_txt(textConnection(input[["text_area"]]))
    })
    
    if(exists("input_sequences")) {
      if(length(input_sequences) > 50) {
        #dummy error, just to stop further processing
        stop("Too many sequences.")
      } else {
        pred_vals(pred_list[[input[["seq_type"]]]], 
                  input_sequences, 
                  names(input_sequences))
      }
    } else {
      NULL
    }
  })
  
  output$dynamic_ui <- renderUI({
    if (!is.null(input[["seq_file"]]))
      input_sequences <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({
      if (!is.null(input[["text_area"]]))
        if(input[["text_area"]] != "")
          input_sequences <- read_txt(textConnection(input[["text_area"]]))
    })
    
  })
  
  output$pred_table <- DT::renderDataTable({
    pred_df <- prediction()
    
    formatRound(my_DT(pred_df), 1L:nrow(pred_df), 2)
  })
  
  
  
  output$dynamic_tabset <- renderUI({
    if(is.null(prediction())) {
        tabPanel(title = "",
                 checkboxGroupInput(inputId = "seq_type", "Choose sequence type", 
                                    choices = c(`16S rRNA` = "rna", mcrA = "mcra"),
                                    selected = "rna"),
                 textAreaInput(inputId = "text_area", 
                               label = "",
                               width = "100%",
                               rows = 15,
                               placeholder = "Paste sequences (FASTA format required) here...", 
                               resize = "horizontal"),
                 p(""),
                 actionButton("use_area", "Submit data from the field above"),
                 fileInput('seq_file', 'Submit .fasta or .txt file:'))
 
    } else {
      tabPanel(title = "Sequence output",
               DT::dataTableOutput("pred_table"),
               tags$p(HTML("<h3><A HREF=\"javascript:history.go(0)\">Start a new query</A></h3>"))
      )
    }
  })
  
  
})
