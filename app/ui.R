

library(shiny)

shinyUI(fluidPage(tags$head(includeScript("ga.js")),
  #tags$style(includeCSS("./www/report.css")),
  theme = shinythemes::shinytheme("flatly"),
  tags$style(HTML("                  
                  .shiny-input-container:not(.shiny-input-container-inline) {
                  width: 100%;
                  }
                  
                  pre{
                  background: white;
                  }")),
  title = "MethaGramPredictor",
  
  headerPanel(""),
  
  sidebarLayout(
    sidebarPanel(style = "background-color: #e0e0e0;",
                 includeMarkdown("readme.md"),
                 pre(includeText("rna.txt"))
    ),
    
    mainPanel(
      uiOutput("dynamic_tabset")
    )
  )
))