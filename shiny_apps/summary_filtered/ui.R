library(shiny)

# Shiny app to display covariates and metrics summary
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

fluidPage(
  
  # Application title
  titlePanel("summary for metrics/covariates"),
  
  # Sidebar with controls to select a dataset and specify the
  # number of observations to view
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Choose a dataset:", 
                  choices = c("metrics", "covariates")),
      textInput("studyID","Study ID:",""),
      textInput("siteID","Site ID:", ""),
      numericInput("obs", "Number of observations to view:", 10),
      textInput("cov_filter","Covariate filter (same syntax as in Google Sheets)",""),
      textInput("metr_name","Metric name (necessary to visualize metrics summary):","")
      
    ),
    
    # Show a summary of the dataset and an HTML table with the 
    # requested number of observations
    mainPanel(
      verbatimTextOutput("summary"),
      tableOutput("view")
    )
  )
)