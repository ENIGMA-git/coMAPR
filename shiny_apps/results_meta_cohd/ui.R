library(shiny)

# Shiny app to export results of Cohen's D meta-analysis
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Meta Cohens'D analysis results"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("studyID","Study ID",""),
      textInput("lm", "Linear Model", ""),
      textInput("metric", "Metric", "")
      ),

    mainPanel(
      tableOutput("meta_cohd_summary")
    )
  )
))