library(shiny)
library(datasets)
library(RMySQL)
library(dplyr)
library(tidyr)
library(stringr)

# Shiny app to display covariates and metrics summary
# Dmitry Isaev
# Imaging Genetics Center, Institute of Neuroimaging and Informatics (INI), 
# University of Southern California
# 2018



function(input, output) {
  # !!! replace this path to path to your 'lib_eregr_core.R'
  suppressPackageStartupMessages(source(file='/path_to_lib_eregr_core/lib_eregr_core.R'))
  # !!! replace this path to path to your 'database_connect.R'
  suppressPackageStartupMessages(source(file='/path_to_database_connect/database_connect.R'))
  
  prepare_covars <- function() {
    
    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
      query <- sprintf("SELECT * FROM covariates_general 
WHERE covariates_general.session_covar_ID=(SELECT MAX(session_covar_ID) AS cur_sess_covar_ID FROM session_covariates,sites_in_study 
WHERE session_covariates.study_site_ID=sites_in_study.study_site_ID 
AND sites_in_study.studyID='%s'
AND siteID='%s');",input$studyID,input$siteID)
      cov <- dbGetQuery(db_conn,query)
      cov[['cov_value']] <- as.numeric(cov[['cov_value']])    
      cov %>% spread(cov_name,cov_value)    
  }
  prepare_covnames <- function() {

    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
    cov <- dbReadTable(db_conn,name="covariates_general")
    cov_names <- unique(cov$cov_name) %>% unlist()
    cov_names 
  }
  prepare_metrics <- function() {

    db_conn <- database_connect()
    on.exit(dbDisconnect(db_conn))
    query <- sprintf("SELECT metrics_data.session_metr_ID as session_metr_ID ,subjID,metric,ROI,vertex,value,siteID,studyID  FROM metrics_data,session_loadMetrics,sites_in_study
WHERE  metrics_data.session_metr_ID=session_loadMetrics.session_metr_ID
                     AND session_loadMetrics.study_site_ID=sites_in_study.study_site_ID
                     AND metrics_data.metric='%s'
                     AND sites_in_study.studyID='%s'
                     AND sites_in_study.siteID='%s'
                     AND metrics_data.session_metr_ID= ( SELECT MAX(session_loadMetrics.session_metr_ID) AS cur_sess_metr_ID FROM metrics_data,session_loadMetrics,sites_in_study
                     WHERE  metrics_data.session_metr_ID=session_loadMetrics.session_metr_ID
                     AND session_loadMetrics.study_site_ID=sites_in_study.study_site_ID
                     AND metrics_data.metric='%s'
                     AND sites_in_study.studyID='%s'
                     AND sites_in_study.siteID='%s' );",input$metr_name,input$studyID,input$siteID,input$metr_name,input$studyID,input$siteID)
    metr <- dbGetQuery(db_conn,query)
    metr %>% group_by(1) %>% spread(ROI,value) %>% select(-1)
  }
  # Return the requested dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "metrics" = prepare_metrics(),
           "covariates" = {
                      dt<- prepare_covars()
                      if(input$cov_filter!="") {
                        filter_full<- gsub(x = input$cov_filter,pattern = "(__)",replacement = "")
                        dt <- dt %>% filter_(str_c(filter_full))
                      }  
                      dt})
  })
  
  # Generate a summary of the dataset
  output$summary <- renderPrint({
    dataset <- datasetInput()
    summary(dataset)
  })
  
  # Show the first "n" observations
  output$view <- renderTable({
    head(datasetInput(), n = input$obs)
  })
  
}
