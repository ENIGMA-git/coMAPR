library(shiny)

# Shiny app - Quality control of amount of patients/controls in each ROI in the same model
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Results - amount of Patients/Controls per all ROI analyses"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("studyID","Study ID", ""),
      textInput("site","Site",""),
      textInput("lm", "Linear Model", "MDDPatVsCont")
      ),
    
  mainPanel(
      tableOutput("demog_summary")
    )
  )
))