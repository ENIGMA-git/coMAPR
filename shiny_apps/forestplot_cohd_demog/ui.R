library(shiny)

# Shiny app to display forestplots for meta-analysis
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Forestplot"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("studyID","Study ID", ""),
      textInput("roi", "ROI", ""),
      textInput("lm", "Linear Model", ""),
      textInput("metric", "Metric", ""),
      textInput("use_meta", "Show meta-analysis results (T/F)","")
      ),

    # Show forest plot
    mainPanel(
      plotOutput("forestPlot"),
      tableOutput("demog_summary")
    )
  )
))