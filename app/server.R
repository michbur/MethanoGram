
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
            filter = "top", rownames = FALSE)

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
        
        pred_vals(pred_list[["rna"]], 
                  input_sequences, 
                  unname(sapply(input_sequences, attr, which = "Annot")),
                  "rna")
      }
    } else {
      NULL
    }
  })
  
  
  output$pred_table <- DT::renderDataTable({
    pred_df <- prediction()[, -4]
    formatRound(my_DT(pred_df), 3L:ncol(pred_df), 2)
  })
  
  output$benchmark_table <- DT::renderDataTable({
    formatRound(my_DT(read.csv("benchmark_res.csv")), 3L:4, 2)
  })
  
  output$benchmark_table_part <- DT::renderDataTable({
    dat <- filter(read.csv("benchmark_res.csv"), 
                  Input.seq == ifelse(seq_type == "rna", "16S rRNA", "mcrA"))
  
    formatRound(my_DT(dat), 3L:4, 2)
  })
  
  output$dynamic_tabset <- renderUI({
    if(is.null(prediction())) {
      tabsetPanel(
        tabPanel(title = "Sequence input",
                 textAreaInput(inputId = "text_area", 
                               label = "",
                               width = "100%",
                               rows = 15,
                               placeholder = "Paste sequences (FASTA format required) here...", 
                               resize = "horizontal"),
                 p(""),
                 actionButton("use_area", "Submit data from the field above"),
                 fileInput('seq_file', 'Submit .fasta or .txt file:')
        ),
        tabPanel(title = "Average error",
                 DT::dataTableOutput("benchmark_table")
                 )
      )
    } else {
      tabsetPanel(
        tabPanel(title = "Results",
                 DT::dataTableOutput("pred_table"),
                 tags$p(HTML("<h3><A HREF=\"javascript:history.go(0)\">Start a new query</A></h3>"))
        ),
        tabPanel(title = "Average error",
                 DT::dataTableOutput("benchmark_table")
        )
      )
    }
  })
})
