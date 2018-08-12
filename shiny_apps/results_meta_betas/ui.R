library(shiny)


# Shiny app to export regression coefficients from meta-analysis stage
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyUI(fluidPage(
  
  titlePanel("Meta betas results"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("studyID","Study ID",""),
      textInput("lm", "Linear Model", ""),
      textInput("metric", "Metric", "")
      ),
    
    mainPanel(
      tableOutput("meta_beta_summary")
    )
  )
))