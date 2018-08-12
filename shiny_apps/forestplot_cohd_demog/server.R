library(shiny)
library(datasets)
library(RMySQL)
library(dplyr)
library(tidyr)

# Shiny app to display forestplots for meta-analysis
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018

shinyServer(function(input, output) {
  # !!! replace this path to path to your 'lib_eregr_core.R'
  suppressPackageStartupMessages(source(file='/path_to_lib_eregr_core/lib_eregr_core.R'))
  # !!! replace this path to path to your 'database_connect.R'
  suppressPackageStartupMessages(source(file='/path_to_database_connect/database_connect.R'))
  
  prepare_forestplot <- function() {

    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
    site_list <- dbGetQuery(db_conn,"SELECT * FROM sites_in_study WHERE studyID='MDDDTI'")
    forest_plot_cohd(db_conn,'MDDDTI',input$lm,input$roi,input$metric,site_list$siteID,use_meta=input$use_meta=="T")
    
  }
  prepare_demog_summ <- function() {

    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
    site_list <- dbGetQuery(db_conn,sprintf("SELECT * FROM sites_in_study WHERE studyID='%s'",input$studyID))
    demog_summary(db_conn,input$studyID,input$lm,input$roi,input$metric,site_list$siteID,use_meta=input$use_meta=="T") %>%
        select (-res_keyID,-vertex)
    
  }  
  output$forestPlot <- renderPlot({
    prepare_forestplot()
    })
  output$demog_summary <- renderTable({
    prepare_demog_summ()
    
  })
})
