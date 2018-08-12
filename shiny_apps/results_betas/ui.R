library(shiny)

# Shiny app to output regression coefficients with standard errors
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Results - betas"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("studyID","Study ID",""),
      textInput("site","Site",""),
      textInput("roi", "ROI", ""),
      textInput("lm", "Linear Model", ""),
      textInput("metric", "Metric", "")
      
      ),
    
    mainPanel(
      tableOutput("betas_summary")
    )
  )
))